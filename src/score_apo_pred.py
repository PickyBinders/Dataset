from collections import defaultdict
from tqdm import tqdm
from create_dataset import read_df, load_cif_data

def get_alignments(aln_file, link_type=None):
    alignments = defaultdict(dict)
    with open(aln_file) as f:
        for i, line in tqdm(enumerate(f)):
            if i == 0:
                columns = line.strip().split()
                continue
            parts = dict(zip(columns, line.strip().split()))
            q_entry, t_entry = parts["query"].split('.')[0].upper(), parts["target"]
            if link_type == "afdb":
                t_entry = t_entry.split('.')[0].upper()
            elif link_type == "apo":
                t_entry = t_entry.split('.cif.gz')[0].upper() + "_" + t_entry.split('_')[-1]
            q_chain = parts["query"].split('_')[1]
            if t_entry not in alignments[q_entry]:
                alignments[q_entry][t_entry] = defaultdict(dict)
            alignments[q_entry][t_entry][q_chain] = parts
    return alignments

def get_scores(df_file, cif_data_dir, aln_file, link_type, df_all_file=None, max_hits=10, protein_fident_threshold=0.99, pocket_fident_threshold=1.0):
    alignments = get_alignments(aln_file, link_type=link_type)
    df = read_df(df_file)
    pockets = load_cif_data(cif_data_dir, names_to_load=["pockets"])["pockets"]
    pocket_residues = {}
    single_pocket_ids = set(df[df["num_chains_in_pocket"] == 1]["single_pocket_ID"])
    for single_pocket_id, pocket in tqdm(pockets.items()):
        if single_pocket_id not in single_pocket_ids:
            continue
        pocket_residues[single_pocket_id] = defaultdict(set)
        for residue in pocket["pocket_residues"]:
            pocket_residues[single_pocket_id][residue['chain']].add(residue['residue_index'])
    entry_to_pockets = df[df["num_chains_in_pocket"] == 1].groupby("PDB_ID").apply(lambda x: x["pocket_ID"].unique()).to_dict()
    if link_type == "apo":
        assert df_all_file is not None
        df_all = read_df(df_all_file)
        chain_to_ligand_type = defaultdict(set)
        for i, row in tqdm(df_all.iterrows()):
            chains = row["single_pocket_chains"]
            if str(chains) == "nan":
                continue
            for chain in chains.split("_"):
                chain_to_ligand_type[f"{row['PDB_ID']}_{chain}"].add((row['Ligand'], row['ligand_type'], row['is_artifact']))

    hits_per_pocket = defaultdict(list)
    for q_entry in tqdm(alignments):
        for q_pocket in entry_to_pockets[q_entry]:
            q_entry, q_biounit, q_chain, q_ligand_chains, q_ligands = q_pocket.split("__")
            q_ligand_chains, q_ligands = q_ligand_chains.split("_"), q_ligands.split("_")
            q_pocket_residues = set()
            for q_ligand_chain in q_ligand_chains:
                q_pocket_residues |= pocket_residues[f"{q_entry}__{q_biounit}__{q_ligand_chain}"][q_chain]
            q_pocket_length = len(q_pocket_residues)
            for target_name in alignments[q_entry]:
                aln = alignments[q_entry][target_name].get(q_chain, None)
                if aln is None:
                    continue
                pocket_lddt = 0
                pocket_fident_query = 0
                q_i, t_i = int(aln['qstart']), int(aln['tstart'])
                for q_a, t_a, lddt in zip(aln['qaln'], aln['taln'], aln['lddtfull'].split(",")):
                    if q_i in q_pocket_residues:
                        if lddt != "nan":
                            pocket_lddt += float(lddt)
                        if t_a == q_a:
                            pocket_fident_query += 1
                    if t_a != "-":
                        t_i += 1
                    if q_a != "-":
                        q_i += 1
                pocket_lddt /= q_pocket_length
                pocket_fident_query /= q_pocket_length
                if pocket_fident_query < pocket_fident_threshold:
                    continue
                if float(aln["fident"]) < protein_fident_threshold:
                    continue
                hit = dict(
                    apo_chain=target_name,
                    protein_lddt_qcov=float(aln["lddt"]) * float(aln["qcov"]),
                    protein_lddt=float(aln["lddt"]),
                    protein_qcov=float(aln["qcov"]),
                    protein_fident=float(aln["fident"]),
                    query_start=int(aln["qstart"]),
                    query_end=int(aln["qend"]),
                    pocket_lddt=pocket_lddt,
                    pocket_fident_query = pocket_fident_query)
                if link_type == "apo":
                    hit["ligands"] = chain_to_ligand_type.get(target_name)
                hits_per_pocket[q_pocket].append(hit)
    hits_per_pocket = {k: sorted(v, key=lambda x: x["protein_qcov"], reverse=True)[:max_hits] for k, v in hits_per_pocket.items()}
    return hits_per_pocket