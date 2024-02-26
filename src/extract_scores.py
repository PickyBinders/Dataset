from tqdm import tqdm
from pathlib import Path
from collections import Counter, defaultdict
import numpy as np
from dataclasses import dataclass
import typing as ty
from create_dataset import read_df, load_cif_data
from extract_plip_data import PLIPHash
import pandas as pnd

PLIP_SUFFIXES = ["", "_nowater"]
INFO_COLUMNS = ["query_pocket", "target_pocket", "protein_mapping", "ligand_mapping"]
SCORE_NAMES = ["protein_lddt", "protein_lddt_qcov", "protein_qcov", "protein_fident", "protein_fident_qcov",
               "pocket_lddt", "pocket_lddt_qcov", "pocket_qcov", "pocket_fident", "pocket_fident_qcov"] + \
                [f"pli_qcov{x}" for x in PLIP_SUFFIXES] + ["ligand_tanimoto_similarity_pocket"]
WEIGHT_TYPES = ["max", "weighted_max", "weighted_sum"]


def combine_scores(q_t_scores):
    q_t_scores_combined = {}
    for s in SCORE_NAMES:
        if s.startswith("ligand"):
            continue
        for suffix in ["_weighted_max", "_weighted_sum"]:
            if f"{s}_foldseek{suffix}" not in q_t_scores and f"{s}_mmseqs{suffix}" not in q_t_scores:
                continue
            s1, s2 = q_t_scores.get(f"{s}_foldseek{suffix}", 0), q_t_scores.get(f"{s}_mmseqs{suffix}", 0)
            source = "foldseek"
            if np.allclose(s1, s2):
                q_t_scores_combined[f"{s}{suffix}"] = s1
                q_t_scores_combined[f"{s}{suffix}_source"] = "both"
            else:
                index = np.argmax([s1, s2])
                q_t_scores_combined[f"{s}{suffix}"] = [s1, s2][index]
                source = ["foldseek", "mmseqs"][index]
                q_t_scores_combined[f"{s}{suffix}_source"] = source
            if not suffix == "_weighted_sum":
                n = f'{s}_{source}{suffix}_mapping'
                if n in q_t_scores:
                    q_t_scores_combined[f"{s}{suffix}_mapping"] = q_t_scores[n]
    return q_t_scores_combined

def load_alignments(source_to_aln_file, chain_mapping, midfix=None, target_has_chain=True):
    data = defaultdict(dict)
    for source, aln_file in source_to_aln_file.items():
        if not aln_file.exists():
            continue
        with open(aln_file) as f:
            for i, line in tqdm(enumerate(f)):
                if i == 0:
                    columns = line.strip().split()
                    continue
                parts = dict(zip(columns, line.strip().split()))
                q_entry, q_chain = parts["query"].replace('.cif.gz', '').replace(".pdb", "").split("_")
                if target_has_chain:
                    t_entry, t_chain = parts["target"].replace('.cif.gz', '').replace(".pdb", "").split("_")
                else:
                    t_entry, t_chain = parts["target"].replace('.cif.gz', '').replace(".pdb", ""), "A"
                if midfix is not None and q_entry[1:3].lower() != midfix:
                    continue
                q_entry, t_entry = q_entry.upper(), t_entry.upper()
                if q_chain not in chain_mapping.get(q_entry, set()) or (target_has_chain and t_chain not in chain_mapping.get(t_entry, set())):
                    continue
                q_chain_mapped = chain_mapping[q_entry][q_chain]
                if target_has_chain:
                    t_chain_mapped = chain_mapping[t_entry][t_chain]
                else:
                    t_chain_mapped = "A"
                if t_entry not in data[q_entry]:
                    data[q_entry][t_entry] = dict()
                if q_chain_mapped not in data[q_entry][t_entry]:
                    data[q_entry][t_entry][q_chain_mapped] = dict()
                if t_chain_mapped not in data[q_entry][t_entry][q_chain_mapped]:
                    data[q_entry][t_entry][q_chain_mapped][t_chain_mapped] = dict()
                lddtfull = [0] * len(parts['qaln'])
                if source == "foldseek":
                    aln_index = 0
                    lddtaln= parts['lddtfull'].split(",")
                    for x, (q_a, t_a) in enumerate(zip(parts['qaln'], parts['taln'])):
                        if q_a != "-" and t_a != "-":
                            lddtfull[x] = float(lddtaln[aln_index])
                            aln_index += 1
                parts['lddtfull'] = lddtfull
                parts["qaln"] = parts["qaln"].upper()
                parts["taln"] = parts["taln"].upper()
                data[q_entry][t_entry][q_chain_mapped][t_chain_mapped][source] = parts
    return data

@dataclass
class PLCScorer:
    entry_to_pockets: ty.Dict[str, ty.List[str]]
    pocket_to_hash: ty.Dict[str, ty.Dict[str, ty.Dict[str, ty.Dict[int, ty.Counter[str]]]]]
    pocket_residues: ty.Dict[str, ty.DefaultDict[str, ty.Dict[int, int]]]
    protein_chain_lengths: ty.Dict[str, ty.Dict[str, int]]
    chain_mapping: ty.Dict[str, ty.Dict[str, str]]
    ligand_to_index: ty.Dict[str, int]
    ligand_to_length: ty.Dict[str, int]
    ligand_similarity_matrix: np.ndarray

    @classmethod
    def initialise(cls, df_file, cif_data_dir, ligand_index_file, ligand_matrix_file, ligand_lengths_file):
        ligand_to_index = pnd.read_csv(ligand_index_file)
        ligand_to_index = dict(zip(ligand_to_index["Ligand"], ligand_to_index.index))
        ligand_to_length = pnd.read_csv(ligand_lengths_file)
        ligand_to_length = dict(zip(ligand_to_length["Ligand"], ligand_to_length["n_ro_bonds"]))
        ligand_similarity_matrix = np.load(ligand_matrix_file)
        cif_data = load_cif_data(cif_data_dir, names_to_load=["pockets", "chain_mapping"])
        df = read_df(df_file)
        df = df[(df["system_ligand_chains_count"] < 10) & (df["system_protein_chains_count"] < 10)].reset_index(drop=True)
        single_pocket_ids = set(df["ligand_system_ID"]).intersection(cif_data["pockets"].keys())
        pocket_residues = {}
        pocket_to_hash = {plip_suffix: {} for plip_suffix in PLIP_SUFFIXES}
        for single_pocket_id, pocket in tqdm(cif_data["pockets"].items()):
            if single_pocket_id not in single_pocket_ids:
                continue
            pocket_residues[single_pocket_id] = defaultdict(dict)
            for plip_suffix in PLIP_SUFFIXES:
                pocket_to_hash[plip_suffix][single_pocket_id] = defaultdict(dict)
            for residue in pocket["pocket_residues"]:
                pocket_residues[single_pocket_id][residue['chain']][residue['residue_index'] + 1] = (residue["residue_number"], residue["one_letter_code"])
        entry_to_pockets = df.groupby("entry_pdb_id").apply(lambda x: x["system_ID"].unique()).to_dict()
        for single_pocket_id, plip_hashes in tqdm(zip(df["ligand_system_ID"], df["ligand_interacting_protein_chains_interactions"])):
            if single_pocket_id not in cif_data["pockets"]:
                continue
            if str(plip_hashes) == "nan":
                continue
            for plip_hash in plip_hashes.split(";"):
                plip_hash = PLIPHash.from_string(plip_hash)
                for plip_suffix in PLIP_SUFFIXES:
                    res = int(plip_hash.residue)
                    if res not in pocket_to_hash[plip_suffix][single_pocket_id][plip_hash.chain]:
                        pocket_to_hash[plip_suffix][single_pocket_id][plip_hash.chain][res] = Counter()
                    h_string = plip_hash.to_string(plip_suffix)
                    if h_string is not None:
                        pocket_to_hash[plip_suffix][single_pocket_id][plip_hash.chain][res] += Counter([h_string])
        protein_chain_lengths = defaultdict(dict)
        for entry, chains, chain_lengths in tqdm(zip(df["entry_pdb_id"], df["ligand_interacting_protein_chains_string"], df["ligand_interacting_protein_chains_lengths"])):
            if str(chains) == "nan":
                continue
            for chain, length in zip(chains.split("_"), chain_lengths):
                if length == "nan":
                    continue
                protein_chain_lengths[entry][chain] = int(length)
        return cls(entry_to_pockets, pocket_to_hash, pocket_residues, protein_chain_lengths, cif_data["chain_mapping"], ligand_to_index, ligand_to_length, ligand_similarity_matrix)
    
    def get_protein_scores_pair(self, aln):
        scores = {}
        for source in aln:
            scores[f"protein_qcov_{source}"] = float(aln[source]["qcov"])
            scores[f"protein_fident_{source}"] = float(aln[source]["fident"])
            scores[f"protein_fident_qcov_{source}"] = float(aln[source]["fident"]) * float(aln[source]["qcov"])
            if source == "foldseek":
                scores[f"protein_lddt_{source}"] = float(aln[source]["lddt"])
                scores[f"protein_lddt_qcov_{source}"] = float(aln[source]["lddt"]) * float(aln[source]["qcov"])
        return scores

    def get_protein_scores(self, q_alignments, q_entry, q_chains, t_entry, t_chains, q_length):
        scores = defaultdict(float)
        mappings = defaultdict(list)
        max_chain_lengths = defaultdict(float)
        protein_chain_mapper = None
        s_matrix = np.zeros((len(q_chains), len(t_chains)))
        for i, q_chain in enumerate(q_chains):
            q_chain_length = self.protein_chain_lengths[q_entry].get(q_chain, 0)
            for j, t_chain in enumerate(t_chains):
                aln = q_alignments.get(t_entry, {}).get(q_chain.split(".")[-1], {}).get(t_chain.split(".")[-1], None)
                if aln is not None:
                    pair_scores = self.get_protein_scores_pair(aln)
                    if protein_chain_mapper is None:
                        if "protein_lddt_qcov_foldseek" in pair_scores:
                            protein_chain_mapper = "protein_lddt_qcov_foldseek"
                        else:
                            protein_chain_mapper = "protein_fident_qcov_mmseqs"
                    s_matrix[i, j] = pair_scores.get(protein_chain_mapper, 0)
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
            aln = q_alignments[t_entry].get(q_chain.split(".")[-1], {}).get(t_chain.split(".")[-1], None)
            if aln is None:
                continue
            q_chain_length = self.protein_chain_lengths[q_entry].get(q_chain, 0)
            pair_scores = self.get_protein_scores_pair(aln)
            for score in pair_scores:
                scores[f"{score}_weighted_sum"] += pair_scores[score] * q_chain_length
            alns[(q_chain, t_chain)] = aln
            mappings[f"{protein_chain_mapper}_weighted_sum"].append((q_chain, t_chain))
            # Zero out the entire row and column for the max element
            s_matrix[q_idx, :] = 0
            s_matrix[:, t_idx] = 0
        for score in scores:
            if score.endswith("_weighted_sum") and q_length > 0:
                scores[score] /= q_length
        return mappings, scores, alns, protein_chain_mapper
    
    def get_pocket_resnum(self, source, i, single_pocket, chain):
        if source == "foldseek":
            return self.pocket_residues[single_pocket][chain].get(i, [None, None])[0]
        else:
            if i in set(x[0] for x in self.pocket_residues[single_pocket][chain].values()):
                return i
        return None

    def get_pocket_pli_scores_pair(self, alns, q_single_pocket, t_single_pocket):
        pocket_scores = defaultdict(float)
        pli_scores = defaultdict(float)
        q_plip = {}
        t_plip = {}
        for q_chain, t_chain in alns:
            for plip_suffix in PLIP_SUFFIXES:
                q_plip[plip_suffix] = self.pocket_to_hash[plip_suffix][q_single_pocket][q_chain]
                t_plip[plip_suffix] = self.pocket_to_hash[plip_suffix][t_single_pocket][t_chain]
            aln = alns[(q_chain, t_chain)]
            for source in aln:
                q_i, t_i = int(aln[source]['qstart']), int(aln[source]['tstart'])
                for i in range(len(aln[source]['qaln'])):
                    q_a, t_a, lddt = aln[source]['qaln'][i], aln[source]['taln'][i], aln[source]['lddtfull'][i]
                    q_n = self.get_pocket_resnum(source, q_i, q_single_pocket, q_chain)
                    if q_n is not None:
                        pocket_scores[f"pocket_lddt_{source}"] += lddt
                        if q_a == t_a:
                            pocket_scores[f"pocket_fident_{source}"] += 1
                        t_n = self.get_pocket_resnum(source, t_i, t_single_pocket, t_chain)
                        if t_n is not None:
                            pocket_scores[f"pocket_qcov_{source}"] += 1
                            pocket_scores[f"pocket_lddt_qcov_{source}"] += lddt
                            if q_a == t_a:
                                pocket_scores[f"pocket_fident_qcov_{source}"] += 1
                            for plip_suffix in PLIP_SUFFIXES:
                                if q_n in q_plip[plip_suffix] and t_n in t_plip[plip_suffix]:
                                    pli_scores[f"pli_qcov{plip_suffix}_{source}"] += sum((q_plip[plip_suffix][q_n] & t_plip[plip_suffix][t_n]).values())
                    if t_a != "-":
                        t_i += 1
                    if q_a != "-":
                        q_i += 1
        return pocket_scores, pli_scores

    def get_pocket_pli_scores(self, alns, q_single_pockets, t_single_pockets, q_pocket_length, q_pli_lengths):
        pocket_scores = defaultdict(float)
        pli_scores = defaultdict(float)
        pocket_mappings = defaultdict(list)
        pli_mappings = defaultdict(list)
        max_pocket_lengths = defaultdict(float)
        max_pli_lengths = defaultdict(float)
        s_matrix = np.zeros((len(q_single_pockets), len(t_single_pockets)))
        ligand_chain_mapper = None

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
                if ligand_chain_mapper is None:
                    if "pocket_lddt_qcov_foldseek" in pocket_pair_scores:
                        ligand_chain_mapper = "pocket_lddt_qcov_foldseek"
                    else:
                        ligand_chain_mapper = "pocket_fident_qcov_mmseqs"
                s_matrix[i, j] = pocket_pair_scores[ligand_chain_mapper] / pocket_length
                for score in pocket_pair_scores:
                    if pocket_pair_scores[score] > pocket_scores[f"{score}_weighted_max"]:
                        pocket_scores[f"{score}_weighted_max"] = pocket_pair_scores[score]
                        max_pocket_lengths[f"{score}_weighted_max"] = pocket_length
                        pocket_mappings[f"{score}_weighted_max"] = [(q_single_pocket, t_single_pocket)]
                    normalized_score = pocket_pair_scores[score] / pocket_length
                    if normalized_score > pocket_scores[f"{score}_max"]:
                        pocket_scores[f"{score}_max"] = normalized_score
                        pocket_mappings[f"{score}_max"] = [(q_single_pocket, t_single_pocket)]
                for score in pli_pair_scores:
                    plip_suffix = score.split("pli_qcov")[-1].split("_")[0]
                    if pli_lengths[plip_suffix] == 0:
                        continue
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
        for score in pli_scores:
            if score.endswith("_weighted_max") and max_pli_lengths[score] > 0:
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
            pocket_mappings[f"{ligand_chain_mapper}_weighted_sum"].append((q_single_pockets[q_idx], t_single_pockets[t_idx]))
            s_matrix[q_idx, :] = 0
            s_matrix[:, t_idx] = 0
        for score in pocket_scores:
            if score.endswith("_weighted_sum"):
                pocket_scores[score] /= q_pocket_length
        for score in pli_scores:
            plip_suffix = score.split("pli_qcov")[-1].split("_")[0]
            if score.endswith("_weighted_sum") and q_pli_lengths[plip_suffix] > 0:
                pli_scores[score] /= q_pli_lengths[plip_suffix]
        return pocket_mappings, pocket_scores, pli_mappings, pli_scores, ligand_chain_mapper

    def get_ligand_scores(self, mapping, q_ligand_length):
        score_name="ligand_tanimoto_similarity_pocket"
        ligand_scores = defaultdict(float)
        ligand_mappings = defaultdict(list)
        max_ligand_lengths = defaultdict(float)

        for (q_ligand, q_ligand_chain), (t_ligand, t_ligand_chain) in mapping:
            if q_ligand not in self.ligand_to_index or t_ligand not in self.ligand_to_index or q_ligand not in self.ligand_to_length or t_ligand not in self.ligand_to_length:
                continue
            ligand_pair_score = self.ligand_similarity_matrix[self.ligand_to_index[q_ligand], self.ligand_to_index[t_ligand]] * self.ligand_to_length[q_ligand]
            if ligand_pair_score > ligand_scores[f"{score_name}_weighted_max"]:
                ligand_scores[f"{score_name}_weighted_max"] = ligand_pair_score
                max_ligand_lengths[f"{score_name}_weighted_max"] = self.ligand_to_length[q_ligand]
                ligand_mappings[f"{score_name}_weighted_max"] = [(q_ligand_chain, t_ligand_chain)]
                normalized_score = ligand_pair_score / self.ligand_to_length[q_ligand]
                if normalized_score > ligand_scores[f"{score_name}_max"]:
                    ligand_scores[f"{score_name}_pocket_max"] = normalized_score
                    ligand_mappings[f"{score_name}_pocket_max"] = [(q_ligand_chain, t_ligand_chain)]
            ligand_scores[f"{score_name}_weighted_sum"] += ligand_pair_score
        ligand_mappings[f"{score_name}_weighted_sum"] = [(q_ligand_chain, t_ligand_chain) for (_, q_ligand_chain), (_, t_ligand_chain) in mapping]
        
        if max_ligand_lengths[f"{score_name}_weighted_max"] > 0:
            ligand_scores[f"{score_name}_weighted_max"] /= max_ligand_lengths[f"{score_name}_weighted_max"]
        if q_ligand_length > 0:
            ligand_scores[f"{score_name}_weighted_sum"] /= q_ligand_length
        return ligand_mappings, ligand_scores

    def get_scores(self, q_entry, q_alignments, suffixes):
        for q_pocket in self.entry_to_pockets[q_entry]:
            q_entry, q_biounit, q_chains, q_ligand_chains, q_ligands = q_pocket.split("__")
            q_chains, q_ligand_chains, q_ligands = q_chains.split("_"), q_ligand_chains.split("_"), q_ligands.split("_")
            q_ligand_chain_to_ligand = dict(zip(q_ligand_chains, q_ligands))
            q_length = sum(self.protein_chain_lengths[q_entry].get(q_chain, 0) for q_chain in q_chains)
            q_single_pockets = [f"{q_entry}__{q_biounit}__{q_ligand_chain}" for q_ligand_chain in q_ligand_chains if f"{q_entry}__{q_biounit}__{q_ligand_chain}" in self.pocket_residues]
            q_pocket_length = sum(len(b) for q_single_pocket in q_single_pockets for b in self.pocket_residues[q_single_pocket].values())
            q_pli_lengths = defaultdict(int)
            for plip_suffix in PLIP_SUFFIXES:
                for q_single_pocket in q_single_pockets:
                    q_pli_lengths[plip_suffix] += sum(sum(self.pocket_to_hash[plip_suffix][q_single_pocket][q_chain][q_res].values()) for q_chain in self.pocket_to_hash[plip_suffix][q_single_pocket] for q_res in self.pocket_to_hash[plip_suffix][q_single_pocket][q_chain])
            q_ligand_length = sum(self.ligand_to_length.get(q_ligand, 0) for q_ligand in q_ligands)
            for t_entry in q_alignments:
                if all(q_chain.split(".")[-1] not in q_alignments[t_entry] for q_chain in q_chains):
                    continue
                all_t_chains = set()
                for q_chain in q_chains:
                    q = q_chain.split(".")[-1]
                    if q in q_alignments[t_entry]:
                        all_t_chains.update(q_alignments[t_entry][q].keys())
                if t_entry not in self.entry_to_pockets:
                    continue
                for t_pocket in self.entry_to_pockets[t_entry]:
                    if t_pocket == q_pocket:
                        continue
                    
                    t_entry, t_biounit, t_chains, t_ligand_chains, t_ligands = t_pocket.split("__")
                    t_chains, t_ligand_chains, t_ligands = t_chains.split("_"), t_ligand_chains.split("_"), t_ligands.split("_")
                    if all(t_chain.split(".")[-1] not in all_t_chains for t_chain in t_chains):
                        continue
                    t_ligand_chain_to_ligand = dict(zip(t_ligand_chains, t_ligands))
                    q_t_scores = {}

                    # Protein score calculation
                    protein_mappings, protein_scores, alns, protein_chain_mapper = self.get_protein_scores(q_alignments, q_entry, q_chains, t_entry, t_chains, q_length)
                    q_t_scores.update(protein_scores)
                    
                    # Pocket and PLI score calculation
                    t_single_pockets = [f"{t_entry}__{t_biounit}__{t_ligand_chain}" for t_ligand_chain in t_ligand_chains if f"{t_entry}__{t_biounit}__{t_ligand_chain}" in self.pocket_residues]
                    pocket_mappings, pocket_scores, pli_mappings, pli_scores, ligand_chain_mapper = self.get_pocket_pli_scores(alns, q_single_pockets, t_single_pockets, q_pocket_length, q_pli_lengths)
                    q_t_scores.update(pocket_scores)
                    q_t_scores.update(pli_scores)
                        
                    # Mapping info
                    for mapping in protein_mappings:
                        q_t_scores[f"{mapping}_mapping"] = ",".join(f"{q_chain}:{t_chain}" for q_chain, t_chain in protein_mappings[mapping])
                    for mapping in pocket_mappings:
                        q_t_scores[f"{mapping}_mapping"] = ",".join(f"{q_ligand_chain.split('__')[-1]}:{t_ligand_chain.split('__')[-1]}" for q_ligand_chain, t_ligand_chain in pocket_mappings[mapping])
                    for mapping in pli_mappings:
                        q_t_scores[f"{mapping}_mapping"] = ",".join(f"{q_ligand_chain.split('__')[-1]}:{t_ligand_chain.split('__')[-1]}" for q_ligand_chain, t_ligand_chain in pli_mappings[mapping])
                    q_t_scores = combine_scores(q_t_scores)
                    if protein_chain_mapper is not None:
                        q_t_scores["protein_mapper"] = "foldseek" if "lddt" in protein_chain_mapper else "mmseqs"
                    if ligand_chain_mapper is not None:
                        q_t_scores["ligand_mapper"] = "foldseek" if "lddt" in ligand_chain_mapper else "mmseqs"
                    q_t_scores["protein_mapping"] = ",".join(f"{q_chain}:{t_chain}" for q_chain, t_chain in protein_mappings[f"{protein_chain_mapper}_weighted_sum"])
                    q_t_scores["ligand_mapping"] = ",".join(f"{q_ligand_chain.split('__')[-1]}:{t_ligand_chain.split('__')[-1]}" for q_ligand_chain, t_ligand_chain in pocket_mappings[f"{ligand_chain_mapper}_weighted_sum"])
                    
                    # Ligand score calculation
                    ligand_mapping_from_pocket = [(q_ligand_chain.split('__')[-1], t_ligand_chain.split('__')[-1]) for q_ligand_chain, t_ligand_chain in pocket_mappings[f"{ligand_chain_mapper}_weighted_sum"]]
                    ligand_mapping_from_pocket = [((q_ligand_chain_to_ligand[q_ligand_chain], q_ligand_chain), (t_ligand_chain_to_ligand[t_ligand_chain], t_ligand_chain)) for q_ligand_chain, t_ligand_chain in ligand_mapping_from_pocket]
                    ligand_mappings, ligand_scores = self.get_ligand_scores(ligand_mapping_from_pocket, q_ligand_length)
                    q_t_scores.update(ligand_scores)
                    for mapping in ligand_mappings:
                        q_t_scores[f"{mapping}_mapping"] = ",".join(f"{q_ligand_chain}:{t_ligand_chain}" for q_ligand_chain, t_ligand_chain in ligand_mappings[mapping])
                    
                    q_t_scores["target_pocket"] = t_pocket
                    q_t_scores["query_pocket"] = q_pocket
                    yield q_t_scores

    def write_scores(self, source_to_aln_file, score_file, midfix=None):
        alignments = load_alignments(source_to_aln_file, self.chain_mapping, midfix=midfix)
        columns = INFO_COLUMNS
        suffixes = ["_weighted_sum", "_weighted_max"]
        for suffix in suffixes:
            for s in SCORE_NAMES:
                columns.append(f"{s}{suffix}")
                if not s.startswith("ligand"):
                    columns.append(f"{s}{suffix}_source")
                if suffix != "_weighted_sum" or s == "ligand_tanimoto_similarity":
                    columns.append(f"{s}{suffix}_mapping")
        columns += ["protein_mapper", "ligand_mapper"]
        with open(score_file, "w") as f:
            f.write("\t".join(columns) + "\n")
            for pdb_id in tqdm(alignments):
                if pdb_id not in self.entry_to_pockets:
                    continue
                for score_dict in self.get_scores(pdb_id, alignments[pdb_id], suffixes):
                    values = [score_dict.get(x, np.nan) for x in columns]
                    values = [f"{x:.3f}" if isinstance(x, float) else x for x in values]
                    f.write("\t".join(str(x) for x in values) + "\n")

    def get_pocket_scores_pair_apo_pred(self, alns, q_single_pocket):
        pocket_scores = defaultdict(float)
        for q_chain, t_chain in alns:
            aln = alns[(q_chain, t_chain)]
            for source in aln:
                q_i, t_i = int(aln[source]['qstart']), int(aln[source]['tstart'])
                for i in range(len(aln[source]['qaln'])):
                    q_a, t_a, lddt = aln[source]['qaln'][i], aln[source]['taln'][i], aln[source]['lddtfull'][i]
                    q_n = self.get_pocket_resnum(source, q_i, q_single_pocket, q_chain)
                    if q_n is not None:
                        pocket_scores[f"pocket_lddt_{source}"] += lddt
                        if q_a == t_a:
                            pocket_scores[f"pocket_fident_{source}"] += 1
                    if t_a != "-":
                        t_i += 1
                    if q_a != "-":
                        q_i += 1
        return pocket_scores

    def get_scores_apo_pred(self, q_entry, q_alignments):
        for q_pocket in self.entry_to_pockets[q_entry]:
            q_entry, q_biounit, q_chain, q_ligand_chains, q_ligands = q_pocket.split("__")
            q_chain, q_ligand_chains, q_ligands= q_chain.split("_"), q_ligand_chains.split("_"), q_ligands.split("_")
            if len(q_chain) > 1:
                # TODO: handle this case?
                continue
            q_chain = q_chain[0]
            q_length = self.protein_chain_lengths[q_entry].get(q_chain, 0)
            q_single_pockets = [f"{q_entry}__{q_biounit}__{q_ligand_chain}" for q_ligand_chain in q_ligand_chains if f"{q_entry}__{q_biounit}__{q_ligand_chain}" in self.pocket_residues]
            for t_entry in q_alignments:
                if q_chain.split(".")[1] in q_alignments[t_entry]:
                    for t_chain in q_alignments[t_entry][q_chain.split(".")[1]]:                      
                        q_t_scores = {}

                        # Protein score calculation
                        protein_mappings, protein_scores, alns, protein_chain_mapper = self.get_protein_scores(q_alignments, q_entry, [q_chain], t_entry, [t_chain], q_length)
                        if len(alns) == 0:
                            continue
                        for key in protein_scores:
                            q_t_scores[key] = protein_scores[key]
                        q_t_scores["protein_mapper"] = "foldseek" if "lddt" in protein_chain_mapper else "mmseqs"
                        # Pocket score calculation
                        pocket_scores = defaultdict(list)
                        for q_single_pocket in q_single_pockets:
                            if q_single_pocket not in self.pocket_residues:
                                continue
                            pocket_length = sum(len(p) for p in self.pocket_residues[q_single_pocket].values())
                            pocket_pair_scores = self.get_pocket_scores_pair_apo_pred(alns, q_single_pocket)
                            for k in pocket_pair_scores:
                                pocket_scores[k].append((q_single_pocket, pocket_pair_scores[k], pocket_length))
                        for k in pocket_scores:
                            sv = sorted(pocket_scores[k], key=lambda x: x[1]/x[2], reverse=True)
                            q_t_scores[f"{k}_weighted_sum"] = sum(x[2] * (x[1] / x[2]) for x in sv) / sum(x[2] for x in sv)
                            q_t_scores[f"{k}_weighted_max"] = sv[0][1] / sv[0][2]
                            q_t_scores[f"{k}_weighted_max_mapping"] = sv[0][0]
                            max_ = sorted(pocket_scores[k], key=lambda x: x[1], reverse=True)[0]
                            q_t_scores[f"{k}_max"] = max_[1] / max_[2]
                            q_t_scores[f"{k}_max_mapping"] = max_[0]
                        for mapping in protein_mappings:
                            q_t_scores[f"{mapping}_mapping"] = ",".join(f"{q_chain}:{t_chain}" for q_chain, t_chain in protein_mappings[mapping])
                        q_t_scores["target"] = f"{t_entry}_{t_chain}"
                        q_t_scores["query_pocket"] = q_pocket
                        yield q_t_scores

    def write_scores_apo_pred(self, source_to_aln_file, score_file, midfix=None, target_has_chain=True):
        alignments = load_alignments(source_to_aln_file, self.chain_mapping, midfix=midfix, target_has_chain=target_has_chain)
        columns = ["query_pocket", "target"]
        for suffix in ["_weighted_sum", "_weighted_max", "_max"]:
            for source in source_to_aln_file.keys():
                for s in SCORE_NAMES:
                    if source == "mmseqs" and "lddt" in s:
                        continue
                    if s.startswith("pocket") and "qcov" in s or s.startswith("pli"):
                        continue
                    columns.append(f"{s}_{source}{suffix}")
                    if suffix != "_weighted_sum":
                        columns.append(f"{s}_{source}{suffix}_mapping")
        columns += ["protein_mapper"]
        with open(score_file, "w") as f:
            f.write("\t".join(columns) + "\n")
            for pdb_id in tqdm(self.entry_to_pockets):
                if pdb_id not in alignments:
                    continue
                for score_dict in self.get_scores_apo_pred(pdb_id, alignments[pdb_id]):
                    values = [score_dict.get(x, np.nan) for x in columns]
                    values = [f"{x:.3f}" if isinstance(x, float) else x for x in values]
                    f.write("\t".join(str(x) for x in values) + "\n")


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--df_file", type=Path, required=True)
    parser.add_argument("--cif_data_dir", type=Path, required=True)
    parser.add_argument("--score_file", type=Path, required=True)
    parser.add_argument("--ligand_index_file", type=Path, required=True)
    parser.add_argument("--ligand_matrix_file", type=Path, required=True)
    parser.add_argument("--ligand_lengths_file", type=Path, required=True)
    parser.add_argument("--aln_file_foldseek", type=Path, default=None)
    parser.add_argument("--aln_file_mmseqs", type=Path, default=None)
    parser.add_argument("--midfix", type=str, default=None)
    parser.add_argument("--target_has_chain", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--apo", action="store_true")
    args = parser.parse_args()
    if not args.score_file.exists() or args.overwrite:
        with open(args.score_file, "w") as f:
            f.write("")
        assert args.aln_file_foldseek is not None or args.aln_file_mmseqs is not None, "one of foldseek or mmseqs aln file required"
        source_to_aln_file = {}
        if args.aln_file_foldseek is not None:
            source_to_aln_file["foldseek"] = args.aln_file_foldseek
        if args.aln_file_mmseqs is not None:
            source_to_aln_file["mmseqs"] = args.aln_file_mmseqs
        if any(x.exists() for x in source_to_aln_file.values()):
            plc = PLCScorer.initialise(args.df_file, args.cif_data_dir, args.ligand_index_file, args.ligand_matrix_file, args.ligand_lengths_file)
            if args.apo:
                plc.write_scores_apo_pred(source_to_aln_file, args.score_file, args.midfix, args.target_has_chain)
            else:
                plc.write_scores(source_to_aln_file, args.score_file, args.midfix)

if __name__ == "__main__":
    main()