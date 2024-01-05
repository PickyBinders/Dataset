from tqdm import tqdm
from pathlib import Path
from collections import Counter, defaultdict
import mmap
import json
import numpy as np
import networkx as nx
from dataclasses import dataclass
import typing as ty
from create_dataset import read_df, load_cif_data

ALN_COLUMNS = "query,target,qlen,lddt,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov,qaln,taln,lddtfull".split(",")

PLIP_SUFFIXES = ["", "_nowater", "_nothreeletter", "_nothreeletterhydrophobic"]
SCORE_COLUMNS = ["query_pocket", "target_pocket", 
                 "chain_mapping", "ligand_mapping",
                        "protein_lddt",
                        "protein_lddt_qcov",
                        "protein_qcov", 
                        "protein_fident", 
                        "pocket_lddt",
                        "pocket_lddt_qcov",
                        "pocket_qcov", 
                        "pocket_fident"] + [f"plip_weighted_jaccard{x}" for x in PLIP_SUFFIXES] + [f"plip_jaccard{x}" for x in PLIP_SUFFIXES]


def get_alignment_offsets(aln_file, output_file, q_entry_prefix=None):
    data = defaultdict(dict)
    with open(aln_file, 'r') as f:
        current_offset = 0  # Keep track of the current byte offset in the file
        for i, line in enumerate(f):
            bytes_in_line = len(line.encode('utf-8'))
            if i == 0:
                columns = line.split()
                current_offset += bytes_in_line
                continue
            parts = dict(zip(columns, line.strip().split()))
            q_entry, t_entry = parts["query"].split('.')[0].upper(), parts["target"].split('.')[0].upper()
            if q_entry_prefix is not None and not q_entry[2:4] == q_entry_prefix.upper():
                current_offset += bytes_in_line
                continue
            q_chain, t_chain = parts["query"].split('_')[1], parts["target"].split('_')[1]
            if t_entry not in data[q_entry]:
                data[q_entry][t_entry] = defaultdict(dict)
            data[q_entry][t_entry][q_chain][t_chain] = current_offset
            current_offset += bytes_in_line
    with open(output_file, "w") as f:
        json.dump(data, f)


def parse_hash(hash_string, water_bridges=True, three_letter=True, keep_three_letter_hydrophobic=False):
    interactions = hash_string.split(';')
    if not three_letter:
        if keep_three_letter_hydrophobic:
            interactions = [":".join(x.split(":")[:1] + x.split(":")[2:]) if not x.startswith("hydrophobic") else x for x in interactions]
        else:
            interactions = [":".join(x.split(":")[:1] + x.split(":")[2:]) for x in interactions]
    if not water_bridges:
        interactions = [x for x in interactions if not x.startswith("water_bridges")]
    return Counter(interactions)

def get_plip_score(hash1, hash2, weighted=True):
    if not weighted:
        hash1, hash2 = set(hash1), set(hash2)

    intersection, union = hash1 & hash2, hash1 | hash2
    if weighted:
        intersection = sum(intersection.values())
        union = sum(union.values())
    else:
        intersection = len(intersection)
        union = len(union)
    if union == 0:
        return np.nan # if both sets are empty, return nan
    return intersection / union

def parse_plip_hash(hash, ptype):
    if ptype == "": return parse_hash(hash)
    elif ptype == "_nowater": return parse_hash(hash, water_bridges=False)
    elif ptype == "_nothreeletter": return parse_hash(hash, three_letter=False, water_bridges=False)
    elif ptype == "_nothreeletterhydrophobic": return parse_hash(hash, three_letter=False, water_bridges=False, keep_three_letter_hydrophobic=True)


@dataclass
class PLCScorer:
    entry_to_pockets: ty.Dict[str, ty.List[str]]
    pocket_to_hash: ty.Dict[str, ty.Dict[str, ty.Counter]]
    pocket_residues: ty.Dict[str, ty.DefaultDict[str, ty.Set[int]]]
    protein_chain_lengths: ty.Dict[str, ty.Dict[str, int]]

    @classmethod
    def initialise(cls, df_file, cif_data_dir):
        cif_data = load_cif_data(cif_data_dir)
        pocket_residues = {}
        for single_pocket_id, pocket in tqdm(cif_data["pockets"].items()):
            pocket_residues[single_pocket_id] = defaultdict(set)
            for residue in pocket["pocket_residues"]:
                pocket_residues[single_pocket_id][residue['chain']].add(residue['residue_index'])
        df = read_df(df_file)
        entry_to_pockets = df.groupby("PDB_ID").apply(lambda x: x["pocket_ID"].unique()).to_dict()
        pocket_to_hash = {}
        for single_pocket_id, plip_hash in tqdm(zip(df["single_pocket_ID"], df["hash"])):
            pocket_to_hash[single_pocket_id] = {plip_suffix: parse_plip_hash(plip_hash, plip_suffix) for plip_suffix in PLIP_SUFFIXES}
        protein_chain_lengths = defaultdict(dict)
        for entry, chains, chain_lengths in tqdm(zip(df["PDB_ID"], df["single_pocket_chains"], df["protein_chain_lengths"])):
            if str(chains) == "nan":
                continue
            for chain, length in zip(chains.split("_"), chain_lengths):
                if length == "nan":
                    continue
                protein_chain_lengths[entry][chain] = int(length)
        return cls(entry_to_pockets, pocket_to_hash, pocket_residues, protein_chain_lengths)

    def get_alignments(self, offsets, mm):
        alignments = {}
        for t_entry in offsets:
            alignments[t_entry] = defaultdict(dict)
            for q_chain in offsets[t_entry]:
                for t_chain, offset in offsets[t_entry][q_chain].items():
                    mm.seek(offset)
                    line = mm.readline().decode('utf-8')
                    parts = dict(zip(ALN_COLUMNS, line.strip().split()))
                    alignments[t_entry][q_chain][t_chain] = parts
        return alignments

    def get_pocket_scores(self, alns, q_single_pocket, t_single_pocket):
        scores = defaultdict(float)
        for q_chain, t_chain in alns:
            aln = alns[(q_chain, t_chain)]
            q_i, t_i = int(aln['qstart']), int(aln['tstart'])
            for q_a, t_a, lddt in zip(aln['qaln'], aln['taln'], aln['lddtfull'].split(",")):
                if q_i in self.pocket_residues[q_single_pocket][q_chain]:
                    if lddt != "nan":
                        scores["pocket_lddt"] += float(lddt)
                    if t_i in self.pocket_residues[t_single_pocket][t_chain]:
                        scores["pocket_qcov"] += 1
                        if lddt != "nan":
                            scores["pocket_lddt_qcov"] += float(lddt)
                        if q_a == t_a:
                            scores["pocket_fident"] += 1
                if t_a != "-":
                    t_i += 1
                if q_a != "-":
                    q_i += 1
        return scores

    @staticmethod
    def greedy_map_chains(q_alignments, q_chains, t_entry, t_chains):
        if len(q_chains) == len(t_chains) == 1:
            return [(q_chains[0], t_chains[0])]
        # Initialize similarity matrix with zeros
        s_matrix = np.zeros((len(q_chains), len(t_chains)))

        # Populate the similarity matrix with lddt * qcov values
        for i, q_chain in enumerate(q_chains):
            for j, t_chain in enumerate(t_chains):
                aln = q_alignments.get(t_entry, {}).get(q_chain, {}).get(t_chain, None)
                if aln is not None:
                    s_matrix[i, j] = float(aln['qcov']) * float(aln['lddt'])

        # do chain mapping
        # Find the highest scores along with their indices in the matrix
        chain_mapping = []
        while np.any(s_matrix):
            score = np.amax(s_matrix)
            if score == 0:
                break
            q_idx, t_idx = np.unravel_index(np.argmax(s_matrix), s_matrix.shape)
            chain_mapping.append((q_chains[q_idx], t_chains[t_idx]))
            # Zero out the entire row and column for the max element
            s_matrix[q_idx, :] = 0
            s_matrix[:, t_idx] = 0
        return chain_mapping

    def greedy_map_ligands(self, alns, q_single_pockets, t_single_pockets):
        if len(q_single_pockets) == len(t_single_pockets) == 1:
            return [(q_single_pockets[0], t_single_pockets[0])]

        # Initialize similarity matrix with zeros
        s_matrix = np.zeros((len(q_single_pockets), len(t_single_pockets)))

        # Populate the similarity matrix with lddt * qcov values
        for i, q_single_pocket in enumerate(q_single_pockets):
            if q_single_pocket not in self.pocket_residues:
                continue
            pocket_length = sum(len(p) for p in self.pocket_residues[q_single_pocket].values())
            for j, t_single_pocket in enumerate(t_single_pockets):
                if t_single_pocket not in self.pocket_residues:
                    continue
                s_matrix[i, j] = self.get_pocket_scores(alns, q_single_pocket, t_single_pocket)["pocket_lddt_qcov"] / pocket_length
        
        # do chain mapping
        # Find the highest scores along with their indices in the matrix
        ligand_mapping = []
        while np.any(s_matrix):
            score = np.amax(s_matrix)
            if score == 0:
                break
            q_idx, t_idx = np.unravel_index(np.argmax(s_matrix), s_matrix.shape)
            ligand_mapping.append((q_single_pockets[q_idx], t_single_pockets[t_idx]))
            # Zero out the entire row and column for the max element
            s_matrix[q_idx, :] = 0
            s_matrix[:, t_idx] = 0
        return ligand_mapping

    def get_protein_scores(self, q_alignments, q_entry, t_entry, chain_mapping):
        score_dict = defaultdict(float)
        alns = {}
        for q_chain, t_chain in chain_mapping:
            aln = q_alignments[t_entry].get(q_chain, {}).get(t_chain, None)
            if aln is None:
                continue
            q_chain_length = self.protein_chain_lengths[q_entry].get(q_chain, 0)
            score_dict["protein_lddt_qcov"] += q_chain_length * (float(aln["lddt"]) * float(aln["qcov"]))
            score_dict["protein_lddt"] += q_chain_length * float(aln["lddt"])
            score_dict["protein_qcov"] += q_chain_length * float(aln["qcov"])
            score_dict["protein_fident"] += q_chain_length * float(aln["fident"])
            alns[(q_chain, t_chain)] = aln
        return score_dict, alns

    def get_scores(self, q_entry, offsets, mm):
        q_alignments = self.get_alignments(offsets, mm)
        scores = defaultdict(list)
        for q_pocket in self.entry_to_pockets[q_entry]:
            q_entry, q_biounit, q_chains, q_ligand_chains, q_ligands = q_pocket.split("__")
            q_chains, q_ligand_chains, q_ligands= q_chains.split("_"), q_ligand_chains.split("_"), q_ligands.split("_")
            q_length = sum(self.protein_chain_lengths[q_entry].get(q_chain, 0) for q_chain in q_chains)
            q_single_pockets = [f"{q_entry}__{q_biounit}__{q_ligand_chain}" for q_ligand_chain in q_ligand_chains]
            q_pocket_length = sum(len(b) for q_single_pocket in q_single_pockets for b in self.pocket_residues[q_single_pocket].values())
            for t_entry in q_alignments:
                if any(q_chain in q_alignments[t_entry] for q_chain in q_chains):
                    if t_entry not in self.entry_to_pockets:
                        continue
                    for t_pocket in self.entry_to_pockets[t_entry]:
                        if t_pocket == q_pocket:
                            continue
                        
                        t_entry, t_biounit, t_chains, t_ligand_chains, t_ligands = t_pocket.split("__")
                        t_chains, t_ligand_chains, t_ligands = t_chains.split("_"), t_ligand_chains.split("_"), t_ligands.split("_")

                        q_t_pocket_scores = {}

                        # Protein score calculation
                        chain_mapping = self.greedy_map_chains(q_alignments, q_chains, t_entry, t_chains)
                        protein_scores, alns = self.get_protein_scores(q_alignments, q_entry, t_entry, chain_mapping)
                        if len(alns) == 0:
                            continue
                        for key in protein_scores:
                            q_t_pocket_scores[key] = protein_scores[key] / q_length

                        # Pocket score calculation
                        t_single_pockets = [f"{t_entry}__{t_biounit}__{t_ligand_chain}" for t_ligand_chain in t_ligand_chains]
                        
                        ligand_mapping = self.greedy_map_ligands(alns, q_single_pockets, t_single_pockets)

                        t_plip = {plip_suffix: Counter() for plip_suffix in PLIP_SUFFIXES}
                        q_plip = {plip_suffix: Counter() for plip_suffix in PLIP_SUFFIXES}
                        pocket_scores = defaultdict(float)
                        for q_single_pocket, t_single_pocket in ligand_mapping:
                            for plip_suffix in PLIP_SUFFIXES:
                                q_plip[plip_suffix] += self.pocket_to_hash[q_single_pocket][plip_suffix]
                                t_plip[plip_suffix] += self.pocket_to_hash[t_single_pocket][plip_suffix]
                            for key, score in self.get_pocket_scores(alns, q_single_pocket, t_single_pocket).items():
                                pocket_scores[key] += score
                        for key in pocket_scores:
                            q_t_pocket_scores[key] = pocket_scores[key] / q_pocket_length
                        
                        for plip_suffix in PLIP_SUFFIXES:
                            q_t_pocket_scores["plip_weighted_jaccard" + plip_suffix] = get_plip_score(q_plip[plip_suffix], t_plip[plip_suffix], weighted=True)
                            q_t_pocket_scores["plip_jaccard" + plip_suffix] = get_plip_score(q_plip[plip_suffix], t_plip[plip_suffix], weighted=False)
                        
                        q_t_pocket_scores["chain_mapping"] = ",".join(f"{q_chain}:{t_chain}" for q_chain, t_chain in chain_mapping)
                        ligand_mapping = [(q.split("_")[-1], t.split("_")[-1]) for q, t in ligand_mapping]
                        q_t_pocket_scores["ligand_mapping"] = ",".join(f"{q_ligand_chain}:{t_ligand_chain}" for q_ligand_chain, t_ligand_chain in ligand_mapping)
                        q_t_pocket_scores["target_pocket"] = t_pocket
                        for score in SCORE_COLUMNS[1:]:
                            if score not in q_t_pocket_scores:
                                q_t_pocket_scores[score] = np.nan
                        if all(np.isnan(x) for x in q_t_pocket_scores.values()):
                            continue
                        scores[q_pocket].append(q_t_pocket_scores)
        return scores


def write_scores(df_file, cif_data_dir, aln_file, offsets_file, score_file, pdb_id_prefix=None, overwrite=False):
    offsets_file = Path(offsets_file)
    if overwrite or not offsets_file.exists():
        get_alignment_offsets(aln_file, offsets_file, pdb_id_prefix)
    with open(offsets_file) as f:
        offsets = json.load(f)
    
    plc = PLCScorer.initialise(df_file, cif_data_dir)
    with open(aln_file, "r+b") as mmap_file:
        mm = mmap.mmap(mmap_file.fileno(), 0)
        with open(score_file, "w") as f:
            f.write("\t".join(SCORE_COLUMNS) + "\n")
            for pdb_id in tqdm(plc.entry_to_pockets):
                if pdb_id not in offsets:
                    continue
                scores = plc.get_scores(pdb_id, offsets[pdb_id], mm)
                for q_pocket in scores:
                    for line in scores[q_pocket]:
                        f.write("\t".join([q_pocket] + [str(line.get(x, np.nan)) for x in SCORE_COLUMNS[1:]]) + "\n")

def get_connected_components(graph, strong=True):
    if strong:
        components = sorted(nx.strongly_connected_components(graph), key=len, reverse=True)
    else:
        components = sorted(nx.connected_components(graph), key=len, reverse=True)
    key_to_component = {}
    for i, c in enumerate(components):
        for k in c:
            key_to_component[k] = i
    return key_to_component


def label_protein_pocket_clusters(score_files, score_name, thresholds, strong=True):
    thresholds = sorted(thresholds)
    graph = nx.DiGraph()
    for score_file in score_files:
        with open(score_file) as f:
            for i, line in tqdm(enumerate(f)):
                if i == 0:
                    continue
                parts = line.strip().split("\t")
                parts = dict(zip(SCORE_COLUMNS, parts))
                score = float(parts[score_name])
                if score < thresholds[0] or np.isnan(score) or str(score) == "nan":
                    continue
                graph.add_edge(parts["query_pocket"], parts["target_pocket"], score=score)
    key_to_component = {thresholds[0]: get_connected_components(graph, strong=strong)}
    for threshold in tqdm(thresholds[1:]):
        graph.remove_edges_from([e for e in graph.edges(data=True) if e[2]["score"] < threshold])
        key_to_component[threshold] = get_connected_components(graph, strong=strong)
    return key_to_component

def label_df(pocket_df, cluster_dir):
    for cluster_file in tqdm(Path(cluster_dir).iterdir()):
        with open(cluster_file) as f:
            key_to_component = json.load(f)
        score_name, strong = cluster_file.stem.split("__")
        for key in key_to_component:
            col = f"{score_name}__{key}__{strong}_component"
            pocket_df[col] = pocket_df["pocket_ID"].apply(lambda x: key_to_component[key].get(x, None))
            max_component = int(pocket_df[col].max()) + 1
            pocket_df.loc[pocket_df[col].isnull(), col] = [int(max_component) + i for i in range(len(pocket_df[pocket_df[col].isnull()]))]
    return pocket_df


# def main():
#     from sys import argv
#     score_file, score_name = argv[1:]
#     name = Path(score_file).stem.split("_scores")[0]
#     score_thresholds = [0.25, 0.5, 0.7, 0.99]
#     key_to_components = label_protein_pocket_clusters(score_file, score_name, score_thresholds)
#     cluster_dir = DATA_FOLDER / "cluster_files"
#     cluster_dir.mkdir(exist_ok=True)
#     with open(cluster_dir / f"{name}__{score_name}.json", "w") as f:
#         json.dump(key_to_components, f)

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--df_file", type=str, required=True)
    parser.add_argument("--cif_data_dir", type=str, required=True)
    parser.add_argument("--aln_file", type=str, required=True)
    parser.add_argument("--prefix", type=str, required=True)
    parser.add_argument("--score_dir", type=str, required=True)
    parser.add_argument("--overwrite", action="store_true")
    
    args = parser.parse_args()
    score_dir = Path(args.score_dir)
    score_dir.mkdir(exist_ok=True)
    score_file = score_dir / f"{args.prefix}_scores.tsv"
    offsets_file = score_dir / f"{args.prefix}_offsets.json"
    if offsets_file.exists():
        write_scores(df_file=args.df_file, 
                    cif_data_dir=args.cif_data_dir, 
                    aln_file=args.aln_file, 
                    offsets_file=offsets_file,
                    score_file=score_file,
                    pdb_id_prefix=args.prefix,
                    overwrite=args.overwrite)


if __name__ == "__main__":
    main()