from tqdm import tqdm
from pathlib import Path
import pandas as pnd
from collections import Counter, defaultdict
import mmap
import json
from multiprocessing import Pool
import numpy as np
import networkx as nx
from picky_binders.annotate.create_dataset import create_dataset_files
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
from rdkit.Chem import rdMolDescriptors as rdmd


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


def get_ligand_distances(mol_dict):
    fp_list = []
    m_list = []
    for m in tqdm(mol_dict):
        try:
            fp_list.append(rdmd.GetMorganFingerprintAsBitVect(mol_dict[m], 3, nBits=2048))
            m_list.append(m)
        except:
            pass
    dists = []
    nfps = len(fp_list)
    for i in tqdm(range(1,nfps)):
        sims = DataStructs.BulkTanimotoSimilarity(fp_list[i],fp_list[:i])
        dists.extend([1-x for x in sims])
    return dists, m_list

def butina_cluster(dists, nfps, cutoff=0.35):
    mol_clusters = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)
    cluster_id_list = [0] * nfps
    for idx,cluster in enumerate(mol_clusters,1):
        for member in cluster:
            cluster_id_list[member] = idx
    return cluster_id_list

def get_alignment_offsets(aln_file):
    data = defaultdict(dict)
    num = 0
    with open(aln_file, 'r') as f:
        current_offset = 0  # Keep track of the current byte offset in the file
        for i, line in tqdm(enumerate(f)):
            bytes_in_line = len(line.encode('utf-8'))
            if i == 0:
                columns = line.split()
                current_offset += bytes_in_line
                continue
            parts = dict(zip(columns, line.strip().split()))
            q_entry, t_entry = parts["query"].split('.')[0].upper(), parts["target"].split('.')[0].upper()
            q_chain, t_chain = parts["query"].split('_')[1], parts["target"].split('_')[1]
            if t_entry not in data[q_entry]:
                data[q_entry][t_entry] = defaultdict(dict)
            data[q_entry][t_entry][q_chain][t_chain] = current_offset
            current_offset += bytes_in_line
            num += 1
    return data


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


@dataclass
class PLC:
    entry_to_pockets: ty.Dict[str, ty.List[str]]
    pocket_to_hash: ty.Dict[str, ty.Dict[str, ty.Counter]]
    pocket_residues: ty.Dict[str, ty.DefaultDict[str, ty.Set[int]]]
    protein_chain_lengths: ty.Dict[str, ty.Dict[str, int]]
    aln_file: ty.Union[str, Path] = None

    @classmethod
    def initialise(cls, df_file, pocket_dir, aln_file, overwrite=False):
        pocket_residues = {}
        for pocket_file in tqdm(Path(pocket_dir).iterdir()):
            with open(pocket_file) as f:
                for entry, pocket in json.load(f).items():
                    pocket_residues[entry] = defaultdict(set)
                    for residue in pocket["binding_site"]:
                        pocket_residues[entry][residue['chain']].add(residue['residue_index'])
        df = read_df(df_file)
        entry_to_pockets = df.groupby("PDB_ID").apply(lambda x: x["pocket_ID"].unique()).to_dict()
        pocket_to_hash = {}
        for single_pocket_id, plip_hash in tqdm(zip(df["joint_pocket_ID"], df["hash"])):
            pocket_to_hash[single_pocket_id] = {plip_suffix: parse_plip_hash(plip_hash, plip_suffix) for plip_suffix in PLIP_SUFFIXES}
        protein_chain_lengths = defaultdict(dict)
        for single_pocket_id, chain_lengths in tqdm(zip(df["single_pocket_ID"], df["protein_chain_lengths"])):
            entry, chains, _, _ = single_pocket_id.split("__")
            entry, biounit = entry.split("_")
            for chain, length in zip(chains.split("_"), chain_lengths):
                if length == "nan":
                    continue
                protein_chain_lengths[entry][chain] = int(length)
        class_object = cls(entry_to_pockets, pocket_to_hash, pocket_residues, protein_chain_lengths, aln_file)
        return class_object

    def get_alignments(self, offsets):
        alignments = {}
        with open(self.aln_file, "r+b") as f:
            mm = mmap.mmap(f.fileno(), 0)
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
                s_matrix[i, j] = self.get_pocket_scores(alns, q_single_pocket, t_single_pocket)["pocket_qcov"] / pocket_length
        
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

    def get_scores(self, q_entry, offsets):
        q_alignments = self.get_alignments(offsets)
        scores = defaultdict(list)
        for q_pocket in self.entry_to_pockets[q_entry]:
            q_entry, q_chains, q_ligand_chains, q_ligands = q_pocket.split("__")
            q_entry, q_biounit = q_entry.split("_")
            q_chains, q_ligand_chains,  q_ligands= q_chains.split("_"), q_ligand_chains.split("_"), q_ligands.split("_")
            q_length = sum(self.protein_chain_lengths[q_entry].get(q_chain, 0) for q_chain in q_chains)
            q_single_pockets = [f"{q_entry}_{q_ligand}_{q_ligand_chain}" for q_ligand_chain, q_ligand in zip(q_ligand_chains, q_ligands)]
            q_pocket_length = sum(len(b) for q_single_pocket in q_single_pockets for b in self.pocket_residues[q_single_pocket].values())
            for t_entry in q_alignments:
                if any(q_chain in q_alignments[t_entry] for q_chain in q_chains):
                    for t_pocket in self.entry_to_pockets[t_entry]:
                        if t_pocket == q_pocket:
                            continue
                        
                        t_entry, t_chains, t_ligand_chains, t_ligands = t_pocket.split("__")
                        t_entry, t_biounit = t_entry.split("_")
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
                        t_single_pockets = [f"{t_entry}_{t_ligand}_{t_ligand_chain}" for t_ligand_chain, t_ligand in zip(t_ligand_chains, t_ligands)]
                        
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


def write_scores_single(args):
    plc, pdb_id, offsets = args
    scores = plc.get_scores(pdb_id, offsets)
    return scores

def write_scores(df_file, pocket_dir, aln_file, offsets_file, score_file, num_threads=1, overwrite=False):
    plc = PLC.initialise(df_file, pocket_dir, aln_file)
    offsets_file = Path(offsets_file)
    if not overwrite and offsets_file.exists():
        with open(offsets_file) as f:
            offsets = json.load(f)
    if overwrite or not offsets_file.exists():
        offsets = get_alignment_offsets(aln_file)
        with open(offsets_file, "w") as f:
            json.dump(offsets, f)
    input_args = [(plc, pdb_id, offsets[pdb_id]) for pdb_id in tqdm(plc.entry_to_pockets) if pdb_id in offsets]
    with open(score_file, "w") as f:
        f.write("\t".join(SCORE_COLUMNS) + "\n")
        with Pool(num_threads) as pool:
            for scores in tqdm(pool.imap(write_scores_single, input_args), total=len(input_args)):
                for q_pocket in scores:
                    for line in scores[q_pocket]:
                        f.write("\t".join([q_pocket] + [str(line.get(x, np.nan)) for x in SCORE_COLUMNS[1:]]) + "\n")



def get_strongly_connected_components(graph):
    components = sorted(nx.strongly_connected_components(graph), key=len, reverse=True)
    key_to_component = {}
    for i, c in enumerate(components):
        for k in c:
            key_to_component[k] = i
    return key_to_component


def label_protein_pocket_clusters(cluster_file, score_file, score_name, thresholds = (0.25, 0.5, 0.7, 0.99)):
    thresholds = sorted(thresholds)
    graph = nx.DiGraph()
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
    key_to_component = {thresholds[0]: get_strongly_connected_components(graph)}
    for threshold in tqdm(thresholds[1:]):
        graph.remove_edges_from([e for e in graph.edges(data=True) if e[2]["score"] < threshold])
        key_to_component[threshold] = get_strongly_connected_components(graph)
    with open(cluster_file, "w") as f:
        json.dump(key_to_component, f)

def label_df(pocket_df, cluster_dir):
    for cluster_file in tqdm(Path(cluster_dir).iterdir()):
        with open(cluster_file) as f:
            key_to_component = json.load(f)
        score_name = cluster_file.stem.split("__")[1]
        for key in key_to_component:
            pocket_df[f"{score_name}__{key}__component"] = pocket_df["pocket_ID"].apply(lambda x: key_to_component[key].get(x, None))
            max_component = int(pocket_df[f"{score_name}__{key}__component"].max()) + 1
            pocket_df.loc[pocket_df[f"{score_name}__{key}__component"].isnull(), f"{score_name}__{key}__component"] = [int(max_component) + i for i in range(len(pocket_df[pocket_df[f"{score_name}__{key}__component"].isnull()]))]
    return pocket_df


# def main():
#     from sys import argv
#     cluster_file, score_file, score_name = argv[1:]
#     score_thresholds = [0.25, 0.5, 0.7, 0.99]
#     key_to_components = label_protein_pocket_clusters(score_file, score_name, score_thresholds)
#     with open(cluster_file, "w") as f:
#         json.dump(key_to_components, f)


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--df_file", type=str, required=True)
    parser.add_argument("--pocket_dir", type=str, required=True)
    parser.add_argument("--aln_file", type=str, required=True)
    parser.add_argument("--offsets_file", type=str, required=True)
    parser.add_argument("--score_file", type=str, required=True)
    parser.add_argument("--num_threads", type=int, default=1)
    
    args = parser.parse_args()
    write_scores(args.df_file, args.pocket_dir, args.aln_file, args.offsets_file, args.score_file, num_threads=args.num_threads)


if __name__ == "__main__":
    main()