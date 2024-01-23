from tqdm import tqdm
from pathlib import Path
from collections import defaultdict
import numpy as np
from dataclasses import dataclass
import typing as ty
from create_dataset import read_df, load_cif_data
from extract_plip_data import PLIPHash

ALN_COLUMNS = "query,target,qlen,lddt,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov,qaln,taln,lddtfull".split(",")

PLIP_SUFFIXES = ["", "_residue", "_nowater"]
INFO_COLUMNS = ["query_pocket", "target_pocket", "protein_mapping", "ligand_mapping"]
SCORE_NAMES = ["protein_lddt", "protein_lddt_qcov", "protein_qcov", "protein_fident", "pocket_lddt", "pocket_lddt_qcov", "pocket_qcov", "pocket_fident"] + [f"pli_qcov{x}" for x in PLIP_SUFFIXES]
PROTEIN_CHAIN_MAPPER = "protein_lddt_qcov"
LIGAND_CHAIN_MAPPER = "pocket_lddt_qcov"

def make_subset_files(df, folder, name):
    """
    Make the subset files for Foldseek.
    Splits PDBs into subsets based on the second two letters of the PDB ID.
    """
    num_chains = 0
    folder = Path(folder)
    subset_folder = folder / "subset_files"
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

@dataclass
class PLCScorer:
    entry_to_pockets: ty.Dict[str, ty.List[str]]
    pocket_to_hash: ty.Dict[str, ty.Dict[str, ty.Dict[str, ty.Dict[int, ty.DefaultDict[str, int]]]]]
    pocket_residues: ty.Dict[str, ty.DefaultDict[str, ty.Dict[int, int]]]
    protein_chain_lengths: ty.Dict[str, ty.Dict[str, int]]

    @classmethod
    def initialise(cls, df_file, cif_data_dir):
        cif_data = load_cif_data(cif_data_dir)
        df = read_df(df_file)
        single_pocket_ids = set(df["single_pocket_ID"])
        pocket_residues = {}
        pocket_to_hash = {plip_suffix: {} for plip_suffix in PLIP_SUFFIXES}
        for single_pocket_id, pocket in tqdm(cif_data["pockets"].items()):
            if single_pocket_id not in single_pocket_ids:
                continue
            pocket_residues[single_pocket_id] = defaultdict(dict)
            for plip_suffix in PLIP_SUFFIXES:
                pocket_to_hash[plip_suffix][single_pocket_id] = defaultdict(dict)
            for residue in pocket["pocket_residues"]:
                pocket_residues[single_pocket_id][residue['chain']][residue['residue_index']] = residue["residue_number"]
        entry_to_pockets = df.groupby("PDB_ID").apply(lambda x: x["pocket_ID"].unique()).to_dict()
        for single_pocket_id, plip_hashes in tqdm(zip(df["single_pocket_ID"], df["hash"])):
            if str(plip_hashes) == "nan":
                continue
            for plip_hash in plip_hashes.split(";"):
                plip_hash = PLIPHash.from_string(plip_hash)
                for plip_suffix in PLIP_SUFFIXES:
                    res = int(plip_hash.residue)
                    if res not in pocket_to_hash[plip_suffix][single_pocket_id][plip_hash.chain]:
                        pocket_to_hash[plip_suffix][single_pocket_id][plip_hash.chain][res] = defaultdict(int)
                    h_string = plip_hash.to_string(plip_suffix)
                    if h_string is not None:
                        pocket_to_hash[plip_suffix][single_pocket_id][plip_hash.chain][res][h_string] += 1
        protein_chain_lengths = defaultdict(dict)
        for entry, chains, chain_lengths in tqdm(zip(df["PDB_ID"], df["single_pocket_chains"], df["protein_chain_lengths"])):
            if str(chains) == "nan":
                continue
            for chain, length in zip(chains.split("_"), chain_lengths):
                if length == "nan":
                    continue
                protein_chain_lengths[entry][chain] = int(length)
        return cls(entry_to_pockets, pocket_to_hash, pocket_residues, protein_chain_lengths)
    
    def get_protein_scores_pair(self, aln):
        scores = {
            "protein_lddt": float(aln["lddt"]),
            "protein_lddt_qcov": float(aln["lddt"]) * float(aln["qcov"]),
            "protein_qcov": float(aln["qcov"]),
            "protein_fident": float(aln["fident"]),
        }
        return scores

    def get_protein_scores(self, q_alignments, q_entry, q_chains, t_entry, t_chains, q_length):
        scores = {}
        mappings = {}
        for score in SCORE_NAMES:
            if score.startswith("protein_"):
                for suffix in ["_weighted_sum", "_weighted_max", "_max"]:
                    scores[score + suffix] = 0
                    mappings[score + suffix] = []
        max_chain_lengths = defaultdict(float)


        s_matrix = np.zeros((len(q_chains), len(t_chains)))
        for i, q_chain in enumerate(q_chains):
            q_chain_length = self.protein_chain_lengths[q_entry].get(q_chain, 0)
            for j, t_chain in enumerate(t_chains):
                aln = q_alignments.get(t_entry, {}).get(q_chain, {}).get(t_chain, None)
                if aln is not None:
                    pair_scores = self.get_protein_scores_pair(aln)
                    s_matrix[i, j] = pair_scores[PROTEIN_CHAIN_MAPPER]
                    for score in pair_scores:
                        if pair_scores[score] * q_chain_length > scores[f"{score}_weighted_max"]:
                            scores[f"{score}_weighted_max"] = pair_scores[score] * q_chain_length
                            max_chain_lengths[f"{score}_weighted_max"] = q_chain_length
                            mappings[f"{score}_weighted_max"] = [(q_chain, t_chain)]
                        if pair_scores[score] > scores[f"{score}_max"]:
                            scores[f"{score}_max"] = pair_scores[score]
                            mappings[f"{score}_max"] = [(q_chain, t_chain)]
        for score in scores:
            if score.endswith("_weighted_max") and max_chain_lengths[score] > 0:
                scores[score] /= max_chain_lengths[score]

        # do chain mapping
        # Find the highest scores along with their indices in the matrix
        alns = {}
        while np.any(s_matrix):
            score = np.amax(s_matrix)
            if score == 0:
                break
            q_idx, t_idx = np.unravel_index(np.argmax(s_matrix), s_matrix.shape)
            q_chain, t_chain = q_chains[q_idx], t_chains[t_idx]
            aln = q_alignments[t_entry].get(q_chain, {}).get(t_chain, None)
            if aln is None:
                continue
            q_chain_length = self.protein_chain_lengths[q_entry].get(q_chain, 0)
            pair_scores = self.get_protein_scores_pair(aln)
            for score in pair_scores:
                scores[f"{score}_weighted_sum"] += pair_scores[score] * q_chain_length
            alns[(q_chain, t_chain)] = aln
            mappings[f"{PROTEIN_CHAIN_MAPPER}_weighted_sum"].append((q_chain, t_chain))
            # Zero out the entire row and column for the max element
            s_matrix[q_idx, :] = 0
            s_matrix[:, t_idx] = 0
        for score in scores:
            if score.endswith("_weighted_sum"):
                scores[score] /= q_length
        return mappings, scores, alns

    def get_pocket_pli_scores_pair(self, alns, q_single_pocket, t_single_pocket):
        pocket_scores = {
            "pocket_lddt": 0,
            "pocket_lddt_qcov": 0,
            "pocket_qcov": 0,
            "pocket_fident": 0,
        }
        pli_scores = {}
        for plip_suffix in PLIP_SUFFIXES:
            pli_scores["pli_qcov" + plip_suffix] = 0
        q_plip = {}
        t_plip = {}
        for q_chain, t_chain in alns:
            for plip_suffix in PLIP_SUFFIXES:
                q_plip[plip_suffix] = self.pocket_to_hash[plip_suffix][q_single_pocket][q_chain]
                t_plip[plip_suffix] = self.pocket_to_hash[plip_suffix][t_single_pocket][t_chain]
            aln = alns[(q_chain, t_chain)]
            q_i, t_i = int(aln['qstart']), int(aln['tstart'])
            for q_a, t_a, lddt in zip(aln['qaln'], aln['taln'], aln['lddtfull'].split(",")):
                q_n = self.pocket_residues[q_single_pocket][q_chain].get(q_i, None)
                if q_n is not None:
                    if lddt != "nan":
                        pocket_scores["pocket_lddt"] += float(lddt)
                    t_n = self.pocket_residues[t_single_pocket][t_chain].get(t_i, None)
                    if t_n is not None:
                        pocket_scores["pocket_qcov"] += 1
                        if lddt != "nan":
                            pocket_scores["pocket_lddt_qcov"] += float(lddt)
                        if q_a == t_a:
                            pocket_scores["pocket_fident"] += 1
                        for plip_suffix in PLIP_SUFFIXES:
                            if q_n in q_plip[plip_suffix] and t_n in t_plip[plip_suffix]:
                                pli_scores[f"pli_qcov{plip_suffix}"] += sum((q_plip[plip_suffix][q_n] & t_plip[plip_suffix][t_n]).values())
                if t_a != "-":
                    t_i += 1
                if q_a != "-":
                    q_i += 1
        return pocket_scores, pli_scores

    def get_pocket_pli_scores(self, alns, q_single_pockets, t_single_pockets, q_pocket_length, q_pli_lengths):
        pocket_scores = {}
        pli_scores = {}
        pocket_mappings = {}
        pli_mappings = {}
        for score in SCORE_NAMES:
            for suffix in ["_weighted_sum", "_weighted_max", "_max"]:
                if score.startswith("pocket_"):
                    pocket_scores[score + suffix] = 0
                    pocket_mappings[score + suffix] = []
                elif score.startswith("pli_"):
                    pli_scores[score + suffix] = 0
                    pli_mappings[score + suffix] = []
        max_pocket_lengths = defaultdict(float)
        max_pli_lengths = defaultdict(float)
        s_matrix = np.zeros((len(q_single_pockets), len(t_single_pockets)))

        for i, q_single_pocket in enumerate(q_single_pockets):
            if q_single_pocket not in self.pocket_residues:
                continue
            pocket_length = sum(len(p) for p in self.pocket_residues[q_single_pocket].values())
            pli_lengths = {}
            for plip_suffix in PLIP_SUFFIXES:
                pli_lengths[plip_suffix] = sum(sum(self.pocket_to_hash[plip_suffix][q_single_pocket][q_chain][q_res].values()) for q_chain in self.pocket_to_hash[plip_suffix][q_single_pocket] for q_res in self.pocket_to_hash[plip_suffix][q_single_pocket][q_chain])
            for j, t_single_pocket in enumerate(t_single_pockets):
                if t_single_pocket not in self.pocket_residues:
                    continue
                pocket_pair_scores, pli_pair_scores = self.get_pocket_pli_scores_pair(alns, q_single_pocket, t_single_pocket)
                s_matrix[i, j] = pocket_pair_scores[LIGAND_CHAIN_MAPPER] / pocket_length
                for score in pocket_pair_scores:
                    if pocket_pair_scores[score] > pocket_scores[f"{score}_weighted_max"]:
                        pocket_scores[f"{score}_weighted_max"] = pocket_pair_scores[score]
                        max_pocket_lengths[f"{score}_weighted_max"] = pocket_length
                        pocket_mappings[f"{score}_weighted_max"] = [(q_single_pocket, t_single_pocket)]
                    if pocket_pair_scores[score] / pocket_length > pocket_scores[f"{score}_max"]:
                        pocket_scores[f"{score}_max"] = pocket_pair_scores[score] / pocket_length
                        pocket_mappings[f"{score}_max"] = [(q_single_pocket, t_single_pocket)]
                for plip_suffix in PLIP_SUFFIXES:
                    score = f"pli_qcov{plip_suffix}"
                    if pli_pair_scores[score] > pli_scores[f"{score}_weighted_max"]:
                        pli_scores[f"{score}_weighted_max"] = pli_pair_scores[score]
                        max_pli_lengths[f"{score}_weighted_max"] = pli_lengths[plip_suffix]
                        pli_mappings[f"{score}_weighted_max"] = [(q_single_pocket, t_single_pocket)]
                    if pli_pair_scores[score] / pli_lengths[plip_suffix] > pli_scores[f"{score}_max"]:
                        pli_scores[f"{score}_max"] = pli_pair_scores[score] / pli_lengths[plip_suffix]
                        pli_mappings[f"{score}_max"] = [(q_single_pocket, t_single_pocket)]
        for score in pocket_scores:
            if score.endswith("_weighted_max") and max_pocket_lengths[score] > 0:
                pocket_scores[score] /= max_pocket_lengths[score]
        for plip_suffix in PLIP_SUFFIXES:
            score = f"pli_qcov{plip_suffix}_weighted_max"
            if max_pli_lengths[score] > 0:
                pli_scores[score] /= max_pli_lengths[score]

        while np.any(s_matrix):
            score = np.amax(s_matrix)
            if score == 0:
                break
            q_idx, t_idx = np.unravel_index(np.argmax(s_matrix), s_matrix.shape)
            q_single_pocket, t_single_pocket = q_single_pockets[q_idx], t_single_pockets[t_idx]
            pocket_pair_scores, pli_pair_scores = self.get_pocket_pli_scores_pair(alns, q_single_pocket, t_single_pocket)
            for score in pocket_pair_scores:
                pocket_scores[f"{score}_weighted_sum"] += pocket_pair_scores[score]
            for score in pli_pair_scores:
                pli_scores[f"{score}_weighted_sum"] += pli_pair_scores[score]
            pocket_mappings[f"{LIGAND_CHAIN_MAPPER}_weighted_sum"].append((q_single_pockets[q_idx], t_single_pockets[t_idx]))
            # Zero out the entire row and column for the max element
            s_matrix[q_idx, :] = 0
            s_matrix[:, t_idx] = 0
        for score in pocket_scores:
            if score.endswith("_weighted_sum"):
                pocket_scores[score] /= q_pocket_length
        for plip_suffix in PLIP_SUFFIXES:
            score = f"pli_qcov{plip_suffix}"
            pli_scores[f"{score}_weighted_sum"] /= q_pli_lengths[plip_suffix]
        return pocket_mappings, pocket_scores, pli_mappings, pli_scores
    

    def get_scores(self, q_entry, q_alignments):
        for q_pocket in self.entry_to_pockets[q_entry]:
            q_entry, q_biounit, q_chains, q_ligand_chains, q_ligands = q_pocket.split("__")
            q_chains, q_ligand_chains, q_ligands= q_chains.split("_"), q_ligand_chains.split("_"), q_ligands.split("_")
            q_length = sum(self.protein_chain_lengths[q_entry].get(q_chain, 0) for q_chain in q_chains)
            q_single_pockets = [f"{q_entry}__{q_biounit}__{q_ligand_chain}" for q_ligand_chain in q_ligand_chains if f"{q_entry}__{q_biounit}__{q_ligand_chain}" in self.pocket_residues]
            q_pocket_length = sum(len(b) for q_single_pocket in q_single_pockets for b in self.pocket_residues[q_single_pocket].values())
            q_pli_lengths = defaultdict(int)
            for plip_suffix in PLIP_SUFFIXES:
                for q_single_pocket in q_single_pockets:
                    q_pli_lengths[plip_suffix] += sum(sum(self.pocket_to_hash[plip_suffix][q_single_pocket][q_chain][q_res].values()) for q_chain in self.pocket_to_hash[plip_suffix][q_single_pocket] for q_res in self.pocket_to_hash[plip_suffix][q_single_pocket][q_chain])
        
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
                        protein_mappings, protein_scores, alns = self.get_protein_scores(q_alignments, q_entry, q_chains, t_entry, t_chains, q_length)
                        if len(alns) == 0:
                            continue
                        for key in protein_scores:
                            q_t_scores[key] = protein_scores[key]

                        # Pocket and PLI score calculation
                        t_single_pockets = [f"{t_entry}__{t_biounit}__{t_ligand_chain}" for t_ligand_chain in t_ligand_chains if f"{t_entry}__{t_biounit}__{t_ligand_chain}" in self.pocket_residues]
                        pocket_mappings, pocket_scores, pli_mappings, pli_scores = self.get_pocket_pli_scores(alns, q_single_pockets, t_single_pockets, q_pocket_length, q_pli_lengths)
                        if len(pocket_mappings[f"{LIGAND_CHAIN_MAPPER}_weighted_sum"]) > 0:
                            for key in pocket_scores:
                                q_t_scores[key] = pocket_scores[key]
                            for key in pli_scores:
                                q_t_scores[key] = pli_scores[key]
                            
                        # Mapping info
                        q_t_scores["protein_mapping"] = ",".join(f"{q_chain}:{t_chain}" for q_chain, t_chain in protein_mappings[f"{PROTEIN_CHAIN_MAPPER}_weighted_sum"])
                        q_t_scores["ligand_mapping"] = ",".join(f"{q_ligand_chain.split('__')[-1]}:{t_ligand_chain.split('__')[-1]}" for q_ligand_chain, t_ligand_chain in pocket_mappings[f"{LIGAND_CHAIN_MAPPER}_weighted_sum"])
                        q_t_scores["target_pocket"] = t_pocket
                        for mapping in protein_mappings:
                            if not mapping.endswith("_weighted_sum"):
                                q_t_scores[f"{mapping}_mapping"] = ",".join(f"{q_chain}:{t_chain}" for q_chain, t_chain in protein_mappings[mapping])
                        for mapping in pocket_mappings:
                            if not mapping.endswith("_weighted_sum"):
                                q_t_scores[f"{mapping}_mapping"] = ",".join(f"{q_ligand_chain.split('__')[-1]}:{t_ligand_chain.split('__')[-1]}" for q_ligand_chain, t_ligand_chain in pocket_mappings[mapping])
                        for mapping in pli_mappings:
                            if not mapping.endswith("_weighted_sum"):
                                q_t_scores[f"{mapping}_mapping"] = ",".join(f"{q_ligand_chain.split('__')[-1]}:{t_ligand_chain.split('__')[-1]}" for q_ligand_chain, t_ligand_chain in pli_mappings[mapping])
                        q_t_scores["query_pocket"] = q_pocket
                        yield q_t_scores
    
    def write_scores(self, aln_file, score_file):
        alignments = load_alignments(aln_file)
        columns = INFO_COLUMNS
        for suffix in ["_weighted_sum", "_weighted_max", "_max"]:
            for s in SCORE_NAMES:
                columns.append(f"{s}{suffix}")
                if suffix != "_weighted_sum":
                    columns.append(f"{s}{suffix}_mapping")
        with open(score_file, "w") as f:
            f.write("\t".join(columns) + "\n")
            for pdb_id in tqdm(self.entry_to_pockets):
                if pdb_id not in alignments:
                    continue
                for score_dict in self.get_scores(pdb_id, alignments[pdb_id]):
                    f.write("\t".join([f"{score_dict.get(x, np.nan)}" if x in INFO_COLUMNS or x.endswith("_mapping") else f"{score_dict.get(x, np.nan):.3f}" for x in columns]) + "\n")


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