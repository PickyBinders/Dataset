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


ALN_COLUMNS = "query,target,lddt,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov,qaln,taln,lddtfull".split(",")

PLIP_SUFFIXES = ["", "_nowater", "_nothreeletter", "_nothreeletterhydrophobic"]
SCORE_COLUMNS = ["query_pocket", "target_pocket", 
                 "chain_mapping", "ligand_mapping",
                        "protein_lddt_qcov", 
                        "protein_lddt",
                        "protein_qcov", 
                        "protein_fident", 
                        "pocket_lddt",
                        "pocket_lddt_shared",
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

def get_chains(rows):
    pdb_id = rows["PDB_ID"].values[0].lower()
    pocket_name = rows["Pocket_Name"].values[0]
    protein_chains = set()
    ligand_chains = set()
    for _, row in rows.iterrows():
        protein_chains |= row["prox_plip_chains"]
        ligand_chains.add(row["ligand_mmcif_chain"])
    return pocket_name, pdb_id, protein_chains, ligand_chains

def get_alignment_offsets(aln_file):
    data = defaultdict(dict)
    num = 0
    with open(aln_file, 'r') as f:
        current_offset = 0  # Keep track of the current byte offset in the file
        pbar = tqdm(enumerate(f))
        for i, line in pbar:
            if i % 100000 == 0:
                pbar.set_description(f"{num} entries")
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


def get_pocket_scores(alns, q_bs, t_bs):
    scores = {"pocket_qcov": 0, "pocket_fident": 0, "pocket_lddt": 0, "pocket_lddt_shared": 0}
    for q_chain, t_chain in alns:
        aln = alns[(q_chain, t_chain)]
        q_i, t_i = int(aln['qstart']), int(aln['tstart'])
        for q_a, t_a, lddt in zip(aln['qaln'], aln['taln'], aln['lddtfull'].split(",")):
            if q_i in q_bs[q_chain]:
                scores["pocket_lddt"] += float(lddt)
            if q_i in q_bs[q_chain] and t_i in t_bs[t_chain]:
                scores["pocket_qcov"] += 1
                scores["pocket_lddt_shared"] += float(lddt)
                if q_a == t_a:
                    scores["pocket_fident"] += 1
            if t_a != "-":
                t_i += 1
            if q_a != "-":
                q_i += 1
    return scores


def greedy_map_chains(q_alignments, q_chains, t_entry, t_chains):
    if len(q_chains) == len(t_chains) == 1:
        return [(q_chains[0], t_chains[0])]

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


def greedy_map_ligands(alns, q_single_pockets, t_single_pockets, binding_sites):
    if len(q_single_pockets) == len(t_single_pockets) == 1:
        return [(q_single_pockets[0], t_single_pockets[0])]

    s_matrix = np.zeros((len(q_single_pockets), len(t_single_pockets)))
    for i, q_single_pocket in enumerate(q_single_pockets):
        if q_single_pocket not in binding_sites:
            continue
        pocket_length = sum(len(b) for b in binding_sites[q_single_pocket].values())
        for j, t_single_pocket in enumerate(t_single_pockets):
            if t_single_pocket not in binding_sites:
                continue
            s_matrix[i, j] = get_pocket_scores(alns, binding_sites[q_single_pocket], binding_sites[t_single_pocket])["pocket_qcov"] / pocket_length
    
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


def get_scores(q_alignments, q_entry, pocket_info, pocket_to_hash, pdb_chain_to_pockets, pdb_chain_lengths):
    scores = defaultdict(list)
    for q_entry, q_biounit, q_chains, q_ligand_chains, q_ligands in pdb_chain_to_pockets[q_entry]:
        q_pocket = f"{q_entry}_{q_biounit}__{'_'.join(q_chains)}__{'_'.join(q_ligand_chains)}__{'_'.join(q_ligands)}"
        for t_entry in q_alignments:
            if any(q_chain in q_alignments[t_entry] for q_chain in q_chains):
                for t_entry, t_biounit, t_chains, t_ligand_chains, t_ligands in pdb_chain_to_pockets[t_entry]:
                    t_pocket = f"{t_entry}_{t_biounit}__{'_'.join(t_chains)}__{'_'.join(t_ligand_chains)}__{'_'.join(t_ligands)}"
                    if t_pocket == q_pocket:
                        continue
                    chain_mapping = greedy_map_chains(q_alignments, q_chains, t_entry, t_chains)
                    score_dict = defaultdict(list)
                    length_dict = defaultdict(list)
                    t_plip = {plip_suffix: Counter() for plip_suffix in PLIP_SUFFIXES}
                    q_plip = {plip_suffix: Counter() for plip_suffix in PLIP_SUFFIXES}
                    alns = {}
                    for q_chain, t_chain in chain_mapping:
                        aln = q_alignments[t_entry].get(q_chain, {}).get(t_chain, None)
                        if aln is None or pdb_chain_lengths[q_entry][q_chain] == "nan" or pdb_chain_lengths[t_entry][t_chain] == "nan":
                            continue
                        score_dict["protein_lddt_qcov"].append(float(aln["lddt"]) * float(aln["qcov"]))
                        score_dict["protein_lddt"].append(float(aln["lddt"]))
                        score_dict["protein_qcov"].append(float(aln["qcov"]))
                        score_dict["protein_fident"].append(float(aln["fident"]))
                        length_dict["protein_length"].append(int(pdb_chain_lengths[q_entry][q_chain]))
                        alns[(q_chain, t_chain)] = aln
                    if not len(alns):
                        continue
                    q_single_pockets = [f"{q_entry}_{q_biounit}__{'_'.join(q_chains)}__{q_ligand_chain}__{q_ligand}" for q_ligand_chain, q_ligand in zip(q_ligand_chains, q_ligands)]
                    t_single_pockets = [f"{t_entry}_{t_biounit}__{'_'.join(t_chains)}__{t_ligand_chain}__{t_ligand}" for t_ligand_chain, t_ligand in zip(t_ligand_chains, t_ligands)]
                    ligand_mapping = greedy_map_ligands(alns, q_single_pockets, t_single_pockets, pocket_info)
                    for q_single_pocket, t_single_pocket in ligand_mapping:
                        length_dict["pocket_length"].append(sum(len(b) for b in pocket_info[q_single_pocket].values()))
                        for plip_suffix in PLIP_SUFFIXES:
                            q_plip[plip_suffix] += pocket_to_hash[q_single_pocket][plip_suffix]
                        for plip_suffix in PLIP_SUFFIXES:
                            t_plip[plip_suffix] += pocket_to_hash[t_single_pocket][plip_suffix]
                        pocket_scores = get_pocket_scores(alns, pocket_info[q_single_pocket], pocket_info[t_single_pocket])
                        for key in pocket_scores:
                            if key.startswith("pocket"):
                                score_dict[key].append(pocket_scores[key])
                    ligand_mapping = [(q.split("__")[2], t.split("__")[2]) for q, t in ligand_mapping]
                    q_t_pocket_scores = {}
                    for key in score_dict:
                        if key.startswith("protein"):
                            try:
                                q_t_pocket_scores[key] = sum(length * score for length, score in zip(length_dict["protein_length"], score_dict[key])) / sum(length_dict["protein_length"])
                            except ZeroDivisionError:
                                q_t_pocket_scores[key] = np.nan
                        elif key.startswith("pocket"):
                            try:
                                q_t_pocket_scores[key] = sum(length * (score / length) for length, score in zip(length_dict["pocket_length"], score_dict[key])) / sum(length_dict["pocket_length"])
                            except ZeroDivisionError:
                                q_t_pocket_scores[key] = np.nan
                    for plip_suffix in PLIP_SUFFIXES:
                        q_t_pocket_scores["plip_weighted_jaccard" + plip_suffix] = get_plip_score(q_plip[plip_suffix], t_plip[plip_suffix], weighted=True)
                        q_t_pocket_scores["plip_jaccard" + plip_suffix] = get_plip_score(q_plip[plip_suffix], t_plip[plip_suffix], weighted=False)
                    q_t_pocket_scores["chain_mapping"] = ",".join(f"{q_chain}:{t_chain}" for q_chain, t_chain in chain_mapping)
                    q_t_pocket_scores["ligand_mapping"] = ",".join(f"{q_ligand_chain}:{t_ligand_chain}" for q_ligand_chain, t_ligand_chain in ligand_mapping)
                    q_t_pocket_scores["target_pocket"] = t_pocket
                    for score in SCORE_COLUMNS[1:]:
                        if score not in q_t_pocket_scores:
                            q_t_pocket_scores[score] = np.nan
                    if all(np.isnan(x) for x in q_t_pocket_scores.values()):
                        continue
                    scores[q_pocket].append(q_t_pocket_scores)
    return scores

def parse_plip_hash(hash, ptype):
    if ptype == "": return parse_hash(hash)
    elif ptype == "_nowater": return parse_hash(hash, water_bridges=False)
    elif ptype == "_nothreeletter": return parse_hash(hash, three_letter=False, water_bridges=False)
    elif ptype == "_nothreeletterhydrophobic": return parse_hash(hash, three_letter=False, water_bridges=False, keep_three_letter_hydrophobic=True)


def get_mappings(df, pocket_dir):
    binding_site_residues = dict()
    for bs_file in tqdm(Path(pocket_dir).iterdir()):
        with open(bs_file) as f:
            binding_site_residues.update(json.load(f))
    pocket_info = {}
    pdb_chain_to_pockets = defaultdict(list)
    pdb_chain_lengths = defaultdict(dict)
    pocket_to_hash = defaultdict(dict)
    pocket_ids = {}
    for _, pocket in tqdm(df.groupby("Pocket_Name")):
        ligand_chains = {}
        protein_chains = []
        for _, row in pocket.iterrows():
            if row["joint_pocket_ID"] not in binding_site_residues:
                continue
            entry, biounit, chains, ligand_chain, ligand = row["PDB_ID"], row["biounit"], sorted(row["prox_plip_chains"]), row["ligand_mmcif_chain"], row["Ligand"]
            biounit = int(biounit)
            for x, chain in enumerate(chains):
                pdb_chain_lengths[entry][chain] = row["protein_chain_lengths"][x]
            single_pocket_id = f"{entry}_{biounit}__{'_'.join(chains)}__{ligand_chain}__{ligand}"
            for p in PLIP_SUFFIXES:
                pocket_to_hash[single_pocket_id][p] = parse_plip_hash(row["hash"], p)
            pocket_info[single_pocket_id] = defaultdict(set)
            for br in binding_site_residues[row["joint_pocket_ID"]]["binding_site"]:
                pocket_info[single_pocket_id][br['chain']].add(br['residue_index'])
            ligand_chains[ligand_chain] = ligand
            protein_chains += chains
        protein_chains = sorted(set(protein_chains))
        ligand_chain_list = sorted(ligand_chains.keys())
        ligands = [ligand_chains[chain] for chain in ligand_chain_list]
        pocket_ids[pocket["Pocket_Name"].values[0]] = f"{entry}_{biounit}__{'_'.join(protein_chains)}__{'_'.join(ligand_chain_list)}__{'_'.join(ligands)}"
        pdb_chain_to_pockets[entry].append((entry, biounit, protein_chains, ligand_chains, ligands))
    return pocket_info, pocket_to_hash, pdb_chain_to_pockets, pdb_chain_lengths, pocket_ids


def write_scores(df_file, aln_file, offsets_file, score_file, num_threads=1):
    df = create_dataset_files.read_df(df_file)
    # keep only first altcode for each pocket
    df = df[df["altcode"].astype(str).isin(["", "A", "nan", " "])].reset_index()
    binding_sites, pocket_to_hash, pdb_chain_to_pockets, pdb_chain_lengths, pocket_ids = get_mappings(df)
    df["pocket_ID"] = df["Pocket_Name"].apply(lambda x: pocket_ids.get(x, None))

    if offsets_file.exists():
        with open(offsets_file) as f:
            offsets = json.load(f)
    if not offsets_file.exists():
        offsets = get_alignment_offsets(aln_file)
        with open(offsets_file, "w") as f:
            json.dump(offsets, f)

    def get_alignments(offsets, q_entry):
        alignments = {}
        if q_entry not in offsets:
            return alignments
        for t_entry in offsets[q_entry]:
            alignments[t_entry] = defaultdict(dict)
            for q_chain in offsets[q_entry][t_entry]:
                for t_chain, offset in offsets[q_entry][t_entry][q_chain].items():
                    MM.seek(offset)
                    line = MM.readline().decode('utf-8')
                    parts = dict(zip(ALN_COLUMNS, line.strip().split()))
                    alignments[t_entry][q_chain][t_chain] = parts
        return alignments

    def write_scores(pdb_id):
        with open(score_file, "a") as f:
            alignments = get_alignments(offsets, pdb_id)
            scores = get_scores(alignments, pdb_id, binding_sites, pocket_to_hash, pdb_chain_to_pockets, pdb_chain_lengths)
            for q_pocket in scores:
                for line in scores[q_pocket]:
                    f.write("\t".join([q_pocket] + [str(line.get(x, np.nan)) for x in SCORE_COLUMNS[1:]]) + "\n")

    with open(aln_file, "r+b") as mmap_f:
        MM = mmap.mmap(mmap_f.fileno(), 0)
        with open(score_file, "w") as f:
            f.write("\t".join(SCORE_COLUMNS) + "\n")

        input_args = [pdb_id for pdb_id in df[df["pocket_ID"].notnull()]["PDB_ID"].unique()]
        with Pool(num_threads) as pool:
            _ = list(tqdm(pool.imap(write_scores, input_args), total=len(input_args)))



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


def main():
    from sys import argv
    cluster_file, score_file, score_name = argv[1:]
    score_thresholds = [0.25, 0.5, 0.7, 0.99]
    key_to_components = label_protein_pocket_clusters(score_file, score_name, score_thresholds)
    with open(cluster_file, "w") as f:
        json.dump(key_to_components, f)


if __name__ == "__main__":
    main()