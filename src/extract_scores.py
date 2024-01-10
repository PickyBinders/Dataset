from tqdm import tqdm
from pathlib import Path
from collections import Counter, defaultdict
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

def make_subset_files(df, folder, name):
    """
    Make the subset files for Foldseek.
    Splits PDBs into subsets based on the second two letters of the PDB ID.
    """
    num_chains = 0
    folder = Path(folder)
    subset_folder = folder / f"{name}_subsets"
    subset_folder.mkdir(exist_ok=True)
    with open(folder / f"{name}.txt" , "w") as f:
        for pdb_id, rows in tqdm(df.groupby("PDB_ID")):
            pdb_id = pdb_id.lower()
            chains = set()
            for _, row in rows.iterrows():
                chains |= set(row["prox_plip_chains"])
            for chain in chains:
                f.write(f"{pdb_id}.cif.gz_{chain}\n")
            prefix = pdb_id[1:3]
            with open(subset_folder / f"{prefix}.txt", "a") as fp:
                for chain in chains:
                    fp.write(f"{pdb_id}.cif.gz_{chain}\n")
            num_chains += len(chains)
    print(f"Number of chains: {num_chains}")


def load_alignments(aln_file):
    data = defaultdict(dict)
    with open(aln_file) as f:
        for i, line in tqdm(enumerate(f)):
            if i == 0:
                columns = line.strip().split()
                continue
            parts = dict(zip(columns, line.strip().split()))
            q_entry, t_entry = parts["query"].split('.')[0].upper(), parts["target"].split('.')[0].upper()
            q_chain, t_chain = parts["query"].split('_')[1], parts["target"].split('_')[1]
            if t_entry not in data[q_entry]:
                data[q_entry][t_entry] = defaultdict(dict)
            data[q_entry][t_entry][q_chain][t_chain] = parts
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
        df = read_df(df_file)
        single_pocket_ids = set(df["single_pocket_ID"])
        pocket_residues = {}
        for single_pocket_id, pocket in tqdm(cif_data["pockets"].items()):
            if single_pocket_id not in single_pocket_ids:
                continue
            pocket_residues[single_pocket_id] = defaultdict(set)
            for residue in pocket["pocket_residues"]:
                pocket_residues[single_pocket_id][residue['chain']].add(residue['residue_index'])
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

    def get_scores(self, q_entry, q_alignments, thresholds=None):
        if thresholds is None:
            thresholds = {}
        for q_pocket in self.entry_to_pockets[q_entry]:
            q_entry, q_biounit, q_chains, q_ligand_chains, q_ligands = q_pocket.split("__")
            q_chains, q_ligand_chains, q_ligands= q_chains.split("_"), q_ligand_chains.split("_"), q_ligands.split("_")
            q_length = sum(self.protein_chain_lengths[q_entry].get(q_chain, 0) for q_chain in q_chains)
            q_single_pockets = [f"{q_entry}__{q_biounit}__{q_ligand_chain}" for q_ligand_chain in q_ligand_chains if f"{q_entry}__{q_biounit}__{q_ligand_chain}" in self.pocket_residues]
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

                        q_t_scores = {}

                        # Protein score calculation
                        chain_mapping = self.greedy_map_chains(q_alignments, q_chains, t_entry, t_chains)
                        protein_scores, alns = self.get_protein_scores(q_alignments, q_entry, t_entry, chain_mapping)
                        if any(protein_scores[x] < thresholds[x] for x in thresholds):
                            continue
                        if len(alns) == 0:
                            continue
                        for key in protein_scores:
                            q_t_scores[key] = protein_scores[key] / q_length

                        # Pocket and PLI score calculation
                        t_single_pockets = [f"{t_entry}__{t_biounit}__{t_ligand_chain}" for t_ligand_chain in t_ligand_chains if f"{t_entry}__{t_biounit}__{t_ligand_chain}" in self.pocket_residues]
                        ligand_mapping = self.greedy_map_ligands(alns, q_single_pockets, t_single_pockets)
                        pocket_scores = defaultdict(float)
                        t_plip = {plip_suffix: Counter() for plip_suffix in PLIP_SUFFIXES}
                        q_plip = {plip_suffix: Counter() for plip_suffix in PLIP_SUFFIXES}
                        for q_single_pocket, t_single_pocket in ligand_mapping:
                            for plip_suffix in PLIP_SUFFIXES:
                                q_plip[plip_suffix] += self.pocket_to_hash[q_single_pocket][plip_suffix]
                                t_plip[plip_suffix] += self.pocket_to_hash[t_single_pocket][plip_suffix]
                            for key, score in self.get_pocket_scores(alns, q_single_pocket, t_single_pocket).items():
                                pocket_scores[key] += score
                        for key in pocket_scores:
                            q_t_scores[key] = pocket_scores[key] / q_pocket_length
                        for plip_suffix in PLIP_SUFFIXES:
                            q_t_scores["plip_weighted_jaccard" + plip_suffix] = get_plip_score(q_plip[plip_suffix], t_plip[plip_suffix], weighted=True)
                            q_t_scores["plip_jaccard" + plip_suffix] = get_plip_score(q_plip[plip_suffix], t_plip[plip_suffix], weighted=False)
                        
                        # Mapping info
                        q_t_scores["chain_mapping"] = ",".join(f"{q_chain}:{t_chain}" for q_chain, t_chain in chain_mapping)
                        q_t_scores["ligand_mapping"] = ",".join(f"{q_ligand_chain.split('__')[-1]}:{t_ligand_chain.split('__')[-1]}" for q_ligand_chain, t_ligand_chain in ligand_mapping)
                        q_t_scores["target_pocket"] = t_pocket
                        for score in SCORE_COLUMNS[1:]:
                            if score not in q_t_scores:
                                q_t_scores[score] = np.nan
                        if all(np.isnan(q_t_scores[x]) for x in SCORE_COLUMNS[4:]):
                            continue
                        q_t_scores["query_pocket"] = q_pocket
                        yield q_t_scores
    
    def write_scores(self, aln_file, score_file):
        alignments = load_alignments(aln_file)
        with open(score_file, "w") as f:
            f.write("\t".join(SCORE_COLUMNS) + "\n")
            for pdb_id in tqdm(self.entry_to_pockets):
                if pdb_id not in alignments:
                    continue
                for score_dict in self.get_scores(pdb_id, alignments[pdb_id]):
                    f.write("\t".join([str(score_dict.get(x, np.nan)) for x in SCORE_COLUMNS[:4]] + [f"{score_dict.get(x, np.nan):.3f}" for x in SCORE_COLUMNS[4:]]) + "\n")


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--df_file", type=Path, required=True)
    parser.add_argument("--cif_data_dir", type=Path, required=True)
    parser.add_argument("--prefix", type=str, required=True)
    parser.add_argument("--aln_dir", type=Path, required=True)
    parser.add_argument("--score_dir", type=str, required=True)
    parser.add_argument("--overwrite", action="store_true")
    
    args = parser.parse_args()
    score_dir = Path(args.score_dir)
    score_dir.mkdir(exist_ok=True)
    score_file = score_dir / f"{args.prefix}_scores.tsv"
    aln_file = args.aln_dir / f"aln_{args.prefix}.tsv"
    if aln_file.exists() and (not score_file.exists() or args.overwrite):
        plc = PLCScorer.initialise(args.df_file, args.cif_data_dir)
        plc.write_scores(aln_file, score_file)


if __name__ == "__main__":
    main()