from multiprocessing import Pool
import pandas as pnd
from tqdm import tqdm
from pathlib import Path
import numpy as np
import pandas as pnd
from collections import defaultdict
import json
from rdkit import Chem
from rdkit import RDLogger
import gzip
from mmcif.io.PdbxReader import PdbxReader

VALIDATION_COLUMNS = [
 'rscc',
 'rsr',
 "entry_r_minus_rfree",
 'entry_rfree',
 'entry_r',
 'entry_clashscore',
 'entry_percent_rama_outliers',
 'entry_mean_b_factor',
 'entry_median_b_factor',
 'ligand_mmcif_chain',
 'Ligand',
 'PDB_chain',
 'PDB_ID',
 "altcode",
 "avgoccu",
 "prox_alternative_configuration_residues_flag",
 'ccd_NumRotatableBonds',
 'perc_prox_rscc',
 'iridium_pass',
 'pass',
 'PocketID',
 ]

COLUMN_RENAME = {
                "avgoccu": "ligand_average_occupancy",
                "ccd_NumRotatableBonds": "ligand_rotatable_bonds_count",
                'hash': 'ligand_interacting_protein_chains_interactions',
                'pass': 'ligand_pass_validation_criteria',
                'perc_prox_rscc': 'ligand_perc_prox_rscc',
                'rscc': 'ligand_rscc',
                'rsr': 'ligand_rsr',
                }

def read_df(df_file, columns=None):
    def parse_list_set(x, use_set=False):
        if str(x) == "nan":
            return None
        else:
            if use_set and x == "set()":
                return set()
            l = x[1:-1].replace("'", "").replace("\\n", "").replace("\n", "").replace(",", "").split(" ")
            if use_set:
                return set(l)
            else:
                return l

    df_file = Path(df_file)
    sep = "\t" if df_file.name.endswith(".tsv") else ","
    df = pnd.read_csv(df_file, sep=sep, low_memory=False, usecols=columns)
    df = df.rename(columns={"Ligand": "ligand_ccd_code", "PDB_ID": "entry_pdb_id"})
    if "ligand_ccd_code" in df.columns:
        df["ligand_ccd_code"] = df["ligand_ccd_code"].apply(lambda x: "SODIUM" if x == "NA" or str(x) == "nan" else x)
    for column in ["ligand_neighboring_ligand_chains", "ligand_interacting_protein_chains", "ligand_interacting_residues", "ligand_neighboring_protein_chains", "ligand_neighboring_residues"]:
        if column in df.columns:
            df[column] = df[column].apply(lambda x: parse_list_set(x, use_set=True))
    for column in ["ligand_center_of_mass", "ligand_interacting_protein_chains_lengths", 
                   "ligand_interacting_protein_chains_uniprot_ids", "num_pdbs_for_uniprots"]:
        if column in df.columns:
            df[column] = df[column].apply(lambda x: parse_list_set(x))
    return df


def parse_cofactors(cofactor_file):
    with open(cofactor_file) as f:
        cofactors_json = json.load(f)
    extra = {
        'Ascorbic acid': ['UU3'],
        'Coenzyme F420': ['6J4', 'F42'],
        'Factor F430': ['F43', 'M43'],
        'Pantetheine': ['PNY'],
        'Pantothenic acids': ['66S','8Q1', 'PAU'],
        'Nicotinamide': ['NCA'],
        'Adenosine nucleotides': ['AMP', 'ATP', "ADP"],
        'Guanosine nucleotides': ['GTP', 'GDP'],
        'Tetrahydrofolic acid': ['MEF'],
        'Lumazine': ['DLZ'],
        'Menaquinone': ['MQ8', 'MQ9', 'MQE'],
        'Heme': ['1CP', 'CP3', 'MMP', 'UP2', 'UP3'],
        'Methanopterin': ['H4M', 'H4Z'],
        'Lipoamide': ['LPM'],
        'Ubiquinone': ['DCQ', 'HQE', 'PLQ'],
        'Pyridoxal': ['PXL', 'UEG'],
        'Siderophores': ['488', 'EB4', 'SE8'],
        'Methanofuran': ['MFN'],
        'Vitamin A': ['BCR', 'ECH', 'EQ3', 'RAW'],
        'Vitamin K1': ['PQN']
    }
    cofactors = set()
    for c in cofactors_json:
        for c_list in cofactors_json[c]:
            cofactors |= set(c_list.get("cofactors", []))
    for c in extra:
        cofactors |= set(extra[c])
    return cofactors

def classify_ligand_rdkit(component):
    """As defined in https://doi.org/10.1021/acs.jcim.3c01573"""
    RDLogger.DisableLog('rdApp.*')
    smarts = {"oligopeptide": Chem.MolFromSmarts("C(=O)C[N;D2,D3]C(=[O;D1])CN"), "oligosaccharide": Chem.MolFromSmarts("O-[C;R0,R1]-[C;R0,R1]-[O,S;R0;D2]-[C;R1]-[O;R1]"), "oligonucleotide": Chem.MolFromSmarts("P(=O)([O-,OH])(OC[C;r5])O[C;r5]")}
    for sm in smarts:
        try:
            if component.mol.HasSubstructMatch(smarts[sm]):
                return sm
        except RuntimeError:
            return "invalid"
    prop_dict = component.physchem_properties
    if not len(prop_dict):
        return "other"
    mw_key = "exactmw"
    if mw_key not in prop_dict:
        mw_key = "amw"
    if prop_dict.get(mw_key, 0) < 300 and prop_dict["CrippenClogP"] < 3 and prop_dict["NumHBD"] <= 3 and prop_dict["NumHBA"] <= 3:
        return "fragment"
    elif prop_dict.get(mw_key, 0) < 500 and prop_dict["CrippenClogP"] < 5 and prop_dict["NumHBD"] <= 5 and prop_dict["NumHBA"] <= 10:
        return "drug-like"
    return "other"

def label_artifacts(df, artifact_file, within_entry_threshold=15, num_prox_plip_residues_threshold=2, common_ligand_threshold=2500):
    """
    From https://doi.org/10.1093/nar/gks966:
    First, if the candidate ligand is in the artifact list and appears >15 times in the same structure file, then it is likely to be crystallization additive and is considered as biologically irrelevant.
    Second, if the number of prox plip residues (i.e. number of contacts) is less than two.

    Additional:
    Third, if the candidate ligand is present over 2,500 times across all pockets and is not a cofactor then it is likely to be a crystallization additive and is considered as biologically irrelevant.
    """
    ignore_ligands = set()
    with open(artifact_file) as f:
        for line in f:
            ignore_ligands.add(line.strip().split("\t")[0])
    df["ligand_in_biolip_artifact_list"] = df["ligand_ccd_code"].apply(lambda x: x in ignore_ligands)
    df["ligand_is_artifact"] = [False] * len(df)
    df.loc[df[df["ligand_ccd_code"].isin(ignore_ligands)].groupby(["entry_pdb_id", "ligand_ccd_code"]).filter(lambda x: len(x) > within_entry_threshold).index, "ligand_is_artifact"] = True
    df["ligand_is_artifact"] = df["ligand_is_artifact"] | (df["ligand_in_biolip_artifact_list"] & (df["ligand_interacting_residues_count"] < num_prox_plip_residues_threshold))
    ligand_counts = dict(zip(df["ligand_ccd_code"].value_counts().index, df["ligand_ccd_code"].value_counts().values))
    common_ligands = set(str(x) for x in ligand_counts if ligand_counts[x] > common_ligand_threshold)
    df["ligand_is_artifact"] = df["ligand_is_artifact"] | ((df["ligand_in_biolip_artifact_list"]) & (df["ligand_ccd_code"].isin(common_ligands)) & (df["ligand_type"] != "cofactor"))
    return df

def label_ligand_types(df, parsed_components, cofactor_file):
    RDLogger.DisableLog('rdApp.*')
    sm_types = dict(zip(df["ligand_ccd_code"], df["ligtype"]))
    cofactors = parse_cofactors(cofactor_file)
    ligand_types = {}
    for lig in df["ligand_ccd_code"].unique():
        lig = str(lig)
        if sm_types.get(lig, None) == "ION":
            ligand_types[lig] = "ion"
        elif lig in cofactors:
            ligand_types[lig] = "cofactor"
        elif lig in parsed_components:
            ligand_types[lig] = classify_ligand_rdkit(parsed_components[lig].component)
        else:
            ligand_types[lig] = "invalid"
    df["ligand_type"] = df["ligand_ccd_code"].apply(lambda x: ligand_types[str(x)])
    return df


def load_cif_data(cif_data_dir, names_to_load=None):
    if names_to_load is None:
        names_to_load = ["dates", "pockets", "uniprot_ids", "chain_mapping", "lengths", "entity_mapping", "label_chain_mapping"]
    cif_data = {x: dict() for x in names_to_load}
    for cif_data_file in tqdm(Path(cif_data_dir).iterdir()):
        with open(cif_data_file) as f:
            d = json.load(f)
            for x in names_to_load:
                cif_data[x].update(d.get(x,{}))
    return cif_data

def get_uniprot_to_chains(chain_to_uniprots):
    uniprot_to_chains = defaultdict(set)
    for chain, uniprot in tqdm(chain_to_uniprots.items()):
        uniprot_to_chains[uniprot].add(chain)
    return uniprot_to_chains

def annotate_chains_and_residues(df, cif_data):
    df["ligand_interacting_protein_chains_count"] = df["ligand_interacting_protein_chains"].apply(lambda x: len(set(x)) if str(x) != "nan" else 0)
    df["ligand_interacting_protein_chains_string"] = df["ligand_interacting_protein_chains"].apply(lambda x: "_".join(sorted(x)) if str(x) != "nan" else np.nan)
    df["ligand_interacting_residues_count"] = df["ligand_interacting_residues"].apply(lambda x: len(set(x)) if str(x) != "nan" else 0)
    df["ligand_neighboring_protein_chains"] = df["ligand_system_ID"].apply(lambda x: set(y["chain"] for y in cif_data["pockets"][x]["pocket_residues"]) if x in cif_data["pockets"] else set())
    df["ligand_neighboring_protein_chains_count"] = df["ligand_neighboring_protein_chains"].apply(lambda x: len(x))
    df["ligand_neighboring_residues"] = df["ligand_system_ID"].apply(lambda x: set(f'{y["chain"]}:{y["residue_number"]}' for y in cif_data["pockets"][x]["pocket_residues"]) if x in cif_data["pockets"] else set())
    df["ligand_center_of_mass"] = df["ligand_system_ID"].apply(lambda x: cif_data["pockets"][x]["center_of_mass"] if x in cif_data["pockets"] else None)
    df["ligand_interacting_protein_chains_interactions_count"] = df["ligand_interacting_protein_chains_interactions"].apply(lambda x: len(str(x).split(";")) if x is not None and str(x) != "nan" else 0)
    df["has_pocket"] = df["ligand_system_ID"].apply(lambda x: x in cif_data["pockets"])
    df["ligand_interacting_protein_chains_uniprot_ids"] = df.apply(lambda row: [cif_data["uniprot_ids"].get(f'{row["entry_pdb_id"]}_{y.split(".")[-1]}', '') for y in row["ligand_interacting_protein_chains_string"].split("_")] if str(row["ligand_interacting_protein_chains_string"]) != "nan" else np.nan, axis=1)
    uniprot_to_chains = get_uniprot_to_chains(cif_data["uniprot_ids"])
    df["num_pdbs_for_uniprots"] = df["ligand_interacting_protein_chains_uniprot_ids"].apply(lambda x: [len(uniprot_to_chains.get(y, [])) for y in x] if str(x) != "nan" else np.nan)
    df["entry_release_date"] = df["entry_pdb_id"].apply(lambda x: cif_data["dates"].get(x.upper(), np.nan))
    df["entry_release_year"] = df["entry_release_date"].apply(lambda x: x.split("-")[0] if str(x) != "nan" else np.nan)
    df["ligand_interacting_protein_chains_lengths"] = df.apply(lambda row: [cif_data["lengths"][row["entry_pdb_id"]].get(y.split(".")[-1], np.nan) for y in row["ligand_interacting_protein_chains_string"].split("_")] if row["entry_pdb_id"] in cif_data["lengths"] and str(row["ligand_interacting_protein_chains_string"]) != "nan" else np.nan, axis=1)
    df["ligand_interacting_protein_chains_entity_ids"] = df.apply(lambda row: [cif_data["entity_mapping"][row["entry_pdb_id"]].get(y.split(".")[-1], np.nan) for y in row["ligand_interacting_protein_chains_string"].split("_")] if row["entry_pdb_id"] in cif_data["entity_mapping"] and str(row["ligand_interacting_protein_chains_string"]) != "nan" else np.nan, axis=1)
    df["ligand_entity_id"] = df.apply(lambda row: cif_data["entity_mapping"][row["entry_pdb_id"]].get(row["ligand_chain"].split(".")[-1], np.nan) if row["entry_pdb_id"] in cif_data["entity_mapping"] else np.nan, axis=1)
    return df


def assign_pocket(df):
    df["system_number"] = [None] * len(df)
    # Assign to same pocket if ligand_neighboring_ligand_chains overlaps with ligand_chain
    for (pdb_id, biounit), group in tqdm(df.groupby(['entry_pdb_id', 'system_biounit'])):
        pocket_number = 0
        pocket_dict = {}
        for index, row in group.sort_values('ligand_neighboring_ligand_chains_count', ascending=False).iterrows():
            prox_ligand_chains = row['ligand_neighboring_ligand_chains']
            if prox_ligand_chains is None:
                prox_ligand_chains = set()
            prox_ligand_chains = set(prox_ligand_chains)
            prox_ligand_chains.add(row['ligand_chain'])
            assigned = False
            for key, value in pocket_dict.items():
                if row['ligand_chain'] in value or len(prox_ligand_chains.intersection(value)) > 0:
                    df.loc[index, 'system_number'] = key
                    pocket_dict[key] |= prox_ligand_chains
                    assigned = True
                    break
            if not assigned:
                pocket_name = f"{pdb_id}_{int(biounit)}_{pocket_number}"
                ligand_chains = prox_ligand_chains.copy()
                ligand_chains.add(row['ligand_chain'])
                pocket_dict[pocket_name] = ligand_chains
                df.loc[index, 'system_number'] = pocket_name
                pocket_number += 1
    df = df[df["system_number"].notnull()].reset_index(drop=True)
    pocket_to_id = {}
    pocket_to_nums = {}
    for _, pocket in tqdm(df.groupby("system_number")):
        ligand_chain_to_name = {}
        for _, row in pocket.iterrows():
            ligand_chain_to_name[row["ligand_chain"]] = row["ligand_ccd_code"]
        ligand_chains = sorted(ligand_chain_to_name)
        ligand_names = [ligand_chain_to_name[x] for x in sorted(ligand_chain_to_name)]
        protein_chains = sorted(set(y for x in pocket["ligand_interacting_protein_chains"] for y in x if str(y) != "nan"))
        pocket_to_id[pocket["system_number"].values[0]] = f"{pocket['entry_pdb_id'].values[0]}__{int(pocket['system_biounit'].values[0])}__{'_'.join(protein_chains)}__{'_'.join(ligand_chains)}__{'_'.join(ligand_names)}"
        pocket_to_nums[pocket["system_number"].values[0]] = (len(ligand_chains), len(protein_chains))
    df["system_ID"] = df["system_number"].apply(lambda x: pocket_to_id[x])
    df["system_ligand_chains_count"] = df["system_number"].apply(lambda x: pocket_to_nums[x][0])
    df["system_protein_chains_count"] = df["system_number"].apply(lambda x: pocket_to_nums[x][1])
    return df

def get_nextgen_mappings(pdb_id, pdb_nextgen_dir):
    per_chain = defaultdict(set)
    cif_file = Path(pdb_nextgen_dir) / pdb_id.lower()[1:3] / f"pdb_0000{pdb_id.lower()}" / f"pdb_0000{pdb_id.lower()}_xyz-enrich.cif.gz"
    entry_info = {}
    covalent_chains = set()
    try:
        with gzip.open(str(cif_file), 'rt', encoding='utf-8') as f:
            data = []
            prd = PdbxReader(f)
            prd.read(data)
        data = data[0]
        entry_info["entry_pdb_id"] = pdb_id
        mappings = [("entry_oligomeric_state", "pdbx_struct_assembly", "oligomeric_details"),
                    ("entry_determination_method", "exptl", "method"),
                    ("entry_keywords", "struct_keywords", "pdbx_keywords"),
                    ("entry_pH", "exptl_crystal_grow", "pH"),
                    ("entry_resolution", "refine", "ls_d_res_high")]
        for key, obj_name, attr_name in mappings:
            x = data.getObj(obj_name)
            if x is not None:
                entry_info[key] = x.getValueOrDefault(attr_name)
        # SIFTS mapping
        xref = data.getObj('pdbx_sifts_xref_db_segments')
        if xref is not None:
            for a in xref:
                if a[2] == "?":
                    continue
                per_chain[(pdb_id, a[1])].add(f'{a[2]}__{a[3]}')
        
        # Label covalent ligands
        struct_conn = data.getObj('struct_conn')
        if struct_conn is not None:
            for a, b in struct_conn.getCombinationCountsWithConditions(['ptnr1_label_asym_id', 'ptnr2_label_asym_id'],
                                                                                      [('conn_type_id', "eq", 'covale')]):
                covalent_chains.add((pdb_id, a))
                covalent_chains.add((pdb_id, b))        
    except Exception:
        pass
    return entry_info, per_chain, covalent_chains

def label_nextgen(df, pdb_nextgen_dir, chain_mapping, num_threads=32):
    input_args = [(pdb_id, pdb_nextgen_dir) for pdb_id in df["entry_pdb_id"].unique()]
    entry_info_all = defaultdict(dict)
    all_sifts_mapping = {}
    all_covalent_mapping = set()
    with Pool(num_threads) as p:
        for entry_info, per_chain, covalent_chains in p.starmap(get_nextgen_mappings, input_args):
            for x in entry_info:
                if x == "entry_pdb_id":
                    continue
                entry_info_all[x][entry_info["entry_pdb_id"]] = entry_info[x]
            all_sifts_mapping.update(per_chain)
            all_covalent_mapping.update(covalent_chains)
    for a in tqdm(entry_info_all):
        df[a] = df["entry_pdb_id"].apply(lambda row: entry_info_all[a].get(row, np.nan))    
    df["ligand_interacting_protein_chains_sifts"] = df.apply(lambda row: [list(all_sifts_mapping.get((row["entry_pdb_id"], chain_mapping[row["entry_pdb_id"]][x.split(".")[-1]]), [])) for x in row["ligand_interacting_protein_chains_string"].split("_")], axis=1)
    df["ligand_has_covalent_bond"] = df.apply(lambda row: (row["entry_pdb_id"], row["ligand_chain"].split(".")[-1]) in all_covalent_mapping, axis=1)
    return df

def make_subset_files(df, foldseek_folder, mmseqs_folder, chain_mapping, name):
    """
    Make the subset files for Foldseek and mmseqs.
    Splits PDBs into subsets based on the second two letters of the PDB ID.
    """
    folders = [foldseek_folder, mmseqs_folder]
    sources = ["foldseek", "mmseqs"]
    for source, folder in zip(sources, folders):
        folder = Path(folder)
        folder.mkdir(exist_ok=True)
        subset_folder = folder / "subset_files"
        subset_folder.mkdir(exist_ok=True)
        first_time = set()
        with open(folder / f"{name}.txt" , "w") as f:
            for pdb_id, rows in tqdm(df.groupby("entry_pdb_id")):
                chains = set()
                for _, row in rows.iterrows():
                    chains |= set(row["ligand_interacting_protein_chains"])
                chains = set(x.split(".")[-1] for x in chains)
                auth_chains = set(chain_mapping[pdb_id][x] for x in chains if x in chain_mapping[pdb_id])
                for chain in auth_chains:
                    if source == "foldseek":
                        f.write(f"{pdb_id.lower()}.cif.gz_{chain}\n")
                    else:
                        f.write(f"{pdb_id.lower()}_{chain}\n")
                prefix = pdb_id[1:3].lower()
                if prefix not in first_time:
                    first_time.add(prefix)
                    open_mode = "w"
                else:
                    open_mode = "a"
                with open(subset_folder / f"{prefix}.txt", open_mode) as fp:
                    for chain in auth_chains:
                        if source == "foldseek":
                            fp.write(f"{pdb_id.lower()}.cif.gz_{chain}\n")
                        else:
                            fp.write(f"{pdb_id.lower()}_{chain}\n")

def create_dataset_files(dataset_dir, foldseek_dir, mmseqs_dir, plip_file, validation_file, cif_data_dir, components_file, cofactor_file, artifact_file, 
                         pdb_nextgen_dir, num_threads=20, max_protein_chains=10, max_ligand_chains=5, overwrite=False):
    from pdbeccdutils.core import ccd_reader
    """
    dataset_dir: directory to save the dataset files
    foldseek_dir: directory to save the Foldseek files
    mmseqs_dir: directory to save the MMSEQS files
    plip_file: file with PLIP results
    validation_file: file with validation results
    cif_data_dir: directory with cif_data files (split by first two characters of PDB ID)
    components_file: file with PDB chemical component dictionary
    cofactor_file: file with cofactors (https://www.ebi.ac.uk/pdbe/api/pdb/compound/cofactors)
    artifact_file: file with artifacts (https://github.com/kad-ecoli/mmCIF2BioLiP/blob/dc9769f286eafc550f799239486ef64450728246/ligand_list)
    pdb_nextgen_dir: directory with PDB NextGen files
    max_protein_chains: maximum number of protein chains to consider
    max_ligand_chains: maximum number of ligand chains to consider
    num_threads: number of threads to use
    overwrite: whether to overwrite existing files
    """
    dataset_dir = Path(dataset_dir)
    dataset_dir.mkdir(exist_ok=True)
    all_pockets_file = dataset_dir / "all_pockets.csv"
    cif_data = load_cif_data(cif_data_dir)
    if overwrite or not all_pockets_file.exists():
        plip_df = read_df(plip_file)
        plip_df["plip_ligand_chain"] = plip_df["ligand_chain"].apply(lambda x: x.split(".")[-1])
        plip_df["joint_pocket_ID"] = plip_df["entry_pdb_id"] + "_" + plip_df["ligand_ccd_code"] + "_" + plip_df["plip_ligand_chain"]
        print("Number of PLIP pockets:", len(plip_df))
        df = read_df(validation_file, columns=VALIDATION_COLUMNS)
        print("Number of validation pockets:", len(df))
        df.rename(columns={"PocketID": "validation_pocket_ID", "ligand_mmcif_chain": "validation_ligand_chain"}, inplace=True)
        df["validation_ligand_chain"] = df["validation_ligand_chain"].apply(lambda x: "NA" if str(x) == "nan" else x)
        df["joint_pocket_ID"] = df["entry_pdb_id"] + "_" + df["ligand_ccd_code"] + "_" + df["validation_ligand_chain"]
        df = df[df["altcode"].isin([" ", "A"])]
        df = pnd.merge(df, plip_df, 
                    on=['joint_pocket_ID', "entry_pdb_id", "ligand_ccd_code"], 
                    how='outer')
        df.rename(columns=COLUMN_RENAME, inplace=True)
        df["ligand_chain"] = df.apply(lambda row: row["ligand_chain"] if str(row["ligand_chain"]) != "nan" else f'1.{row["validation_ligand_chain"]}', axis=1)
        df["system_biounit"] = df["system_biounit"].fillna(1)
        df["ligand_system_ID"] = df.apply(lambda row: f'{row["entry_pdb_id"]}__{int(row["system_biounit"])}__{row["ligand_chain"]}', axis=1)
        print("Annotation of chains and residues")
        df = annotate_chains_and_residues(df, cif_data)
        RDLogger.DisableLog('rdApp.*')
        parsed_components = ccd_reader.read_pdb_components_file(str(components_file))
        print("Labeling ligand types")
        df = label_ligand_types(df, parsed_components, cofactor_file)
        print("Labeling artifacts")
        df = label_artifacts(df, artifact_file)
        df["has_plip"] = df["plip_pocket_ID"].notna()
        df["has_validation"] = df["validation_pocket_ID"].notna()
        df.to_csv(all_pockets_file, index=False)
    df = read_df(all_pockets_file)
    print("Number of pockets after merging PLIP and validation:", len(df), df["ligand_system_ID"].nunique())
    df = df[df["has_plip"] & df["has_pocket"]].reset_index(drop=True)
    df["ligand_neighboring_ligand_chains_count"] = df["ligand_neighboring_ligand_chains"].apply(lambda x: len(x) if x is not None else 0)
    df = assign_pocket(df)
    print("Number of pockets after filtering for PLIP and binding sites and merging pockets")
    print("\tTotal rows:", len(df))
    print("\tTotal pockets:", df["ligand_system_ID"].nunique())
    print("\tTotal merged pockets:", df["system_ID"].nunique())
    df.to_csv(dataset_dir / "all_pockets_with_plip_and_bs.csv", index=False)
    df["ligand_is_ion"] = df["ligand_type"] == "ion"
    only_ions = df.groupby("system_ID").agg({"ligand_is_ion": "all"})
    only_ions = only_ions[only_ions["ligand_is_ion"]].index
    df = df[~df["system_ID"].isin(only_ions)]
    df.to_csv(dataset_dir / "small_molecule_pockets.csv", index=False)
    print("Number of pockets with small molecules:", df["system_ID"].nunique())
    sm_pockets = df[(~df["ligand_is_ion"]) & (~df["ligand_is_artifact"])]["system_ID"].unique()
    df = df[df["system_ID"].isin(sm_pockets)]
    print("Number of pockets with small molecules (no artifacts):", df["system_ID"].nunique())
    df.to_csv(dataset_dir / "small_molecule_pockets_no_artifacts.csv", index=False)
    df = df[(df["system_protein_chains_count"] <= max_protein_chains) & (df["system_ligand_chains_count"] <= max_ligand_chains)]
    print("Number of pockets with small molecules (no artifacts, max protein chains, max ligand chains):", df["system_ID"].nunique())
    print("Annotating info from SMTL")
    df = label_nextgen(df, pdb_nextgen_dir, cif_data["label_chain_mapping"], num_threads=num_threads)
    df.to_csv(dataset_dir / "filtered_pockets.csv", index=False)
    make_subset_files(df, Path(foldseek_dir) / "filtered_pockets", Path(mmseqs_dir) / "filtered_pockets", cif_data["label_chain_mapping"], "filtered_pockets")

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset_dir", type=str, help="Directory to save the dataset files")
    parser.add_argument("foldseek_dir", type=str, help="Directory to save the Foldseek files")
    parser.add_argument("mmseqs_dir", type=str, help="Directory to save the MMSEQS files")
    parser.add_argument("plip_file", type=str, help="File with PLIP results")
    parser.add_argument("validation_file", type=str, help="File with validation results")
    parser.add_argument("cif_data_dir", type=str, help="Directory with cif_data files")
    parser.add_argument("components_file", type=str, help="File with PDB chemical component dictionary")
    parser.add_argument("cofactor_file", type=str, help="File with cofactors (https://www.ebi.ac.uk/pdbe/api/pdb/compound/cofactors)")
    parser.add_argument("artifact_file", type=str, help="File with artifacts (https://github.com/kad-ecoli/mmCIF2BioLiP/blob/dc9769f286eafc550f799239486ef64450728246/ligand_list)")
    parser.add_argument("pdb_nextgen_dir", type=str, help="Directory with PDB NextGen files", default=Path("/scicore/data/managed/PDB_NEXTGEN/latest/pdb_nextgen/data/entries/divided/"))
    parser.add_argument("--num_threads", type=int, default=20, help="Number of threads to use")
    parser.add_argument("--max_protein_chains", type=int, default=5, help="Maximum number of protein chains to consider")
    parser.add_argument("--max_ligand_chains", type=int, default=5, help="Maximum number of ligand chains to consider")
    parser.add_argument("--overwrite", action="store_true", help="Whether to overwrite existing files")

    args = parser.parse_args()
    create_dataset_files(**vars(args))

if __name__ == "__main__":
    main()