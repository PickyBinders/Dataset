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
from pdbeccdutils.core import ccd_reader

VALIDATION_COLUMNS = [
 'rscc',
 'rsr',
 "entry_resolution",
 "entry_r_minus_rfree",
 'entry_rfree',
 'entry_r',
 'entry_clashscore',
 'entry_percent_rama_outliers',
 'entry_mean_b_factor',
 'entry_median_b_factor',
 'ligand_mmcif_chain',
 'ligand_num_clashes',
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
    if "ligand_mmcif_chain" in df.columns:
        df["ligand_mmcif_chain"] = df["ligand_mmcif_chain"].apply(lambda x: "NA" if str(x) == "nan" else x)
    df["Ligand"] = df["Ligand"].apply(lambda x: "SODIUM" if x == "NA" or str(x) == "nan" else x)
    if "ligand_mmcif_chain" in df.columns:
        df["joint_pocket_ID"] = df["PDB_ID"] + "_" + df["Ligand"] + "_" + df["ligand_mmcif_chain"]
    for column in ["prox_plip_ligand_chains", "prox_plip_chains", "prox_plip_residues", "prox_chains", "prox_residues"]:
        if column in df.columns:
            df[column] = df[column].apply(lambda x: parse_list_set(x, use_set=True))
    for column in ["center_of_mass", "protein_chain_lengths", "UniProt_IDs", "num_uniprot_ids", "num_pdbs_for_uniprots"]:
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
    df["is_biolip_artifact"] = df["Ligand"].apply(lambda x: x in ignore_ligands)
    df["is_artifact"] = [False] * len(df)
    df.loc[df[df["Ligand"].isin(ignore_ligands)].groupby(["PDB_ID", "Ligand"]).filter(lambda x: len(x) > within_entry_threshold).index, "is_artifact"] = True
    df["is_artifact"] = df["is_artifact"] | (df["is_biolip_artifact"] & (df["num_prox_plip_residues"] < num_prox_plip_residues_threshold))
    ligand_counts = dict(zip(df["Ligand"].value_counts().index, df["Ligand"].value_counts().values))
    common_ligands = set(str(x) for x in ligand_counts if ligand_counts[x] > common_ligand_threshold)
    df["is_artifact"] = df["is_artifact"] | ((df["is_biolip_artifact"]) & (df["Ligand"].isin(common_ligands)) & (df["ligand_type"] != "cofactor"))
    return df

def label_ligand_types(df, parsed_components, cofactor_file):
    RDLogger.DisableLog('rdApp.*')
    sm_types = dict(zip(df["Ligand"], df["ligtype"]))
    cofactors = parse_cofactors(cofactor_file)
    ligand_types = {}
    for lig in df["Ligand"].unique():
        lig = str(lig)
        if sm_types.get(lig, None) == "ION":
            ligand_types[lig] = "ion"
        elif lig in cofactors:
            ligand_types[lig] = "cofactor"
        elif lig in parsed_components:
            ligand_types[lig] = classify_ligand_rdkit(parsed_components[lig].component)
        else:
            ligand_types[lig] = "invalid"
    df["ligand_type"] = df["Ligand"].apply(lambda x: ligand_types[str(x)])
    return df


def load_cif_data(cif_data_dir, names_to_load=None):
    if names_to_load is None:
        names_to_load = ["dates", "pockets", "uniprot_ids", "chain_mapping", "lengths", "entity_mapping"]
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

def annotate_chains_and_residues(df, cif_data_dir):
    cif_data = load_cif_data(cif_data_dir)
    df["num_prox_plip_chains"] = df["prox_plip_chains"].apply(lambda x: len(set(x)) if str(x) != "nan" else 0)
    df["num_prox_plip_residues"] = df["prox_plip_residues"].apply(lambda x: len(set(x)) if str(x) != "nan" else 0)
    df["prox_chains"] = df["single_pocket_ID"].apply(lambda x: set(y["chain"] for y in cif_data["pockets"][x]["pocket_residues"]) if x in cif_data["pockets"] else set())
    df["num_prox_chains"] = df["prox_chains"].apply(lambda x: len(x))
    df["prox_residues"] = df["single_pocket_ID"].apply(lambda x: set(f'{y["chain"]}:{y["residue_number"]}' for y in cif_data["pockets"][x]["pocket_residues"]) if x in cif_data["pockets"] else set())
    df["center_of_mass"] = df["single_pocket_ID"].apply(lambda x: cif_data["pockets"][x]["center_of_mass"] if x in cif_data["pockets"] else None)
    df["num_plip_interactions"] = df["hash"].apply(lambda x: len(str(x).split(";")) if x is not None and str(x) != "nan" else 0)
    df["num_unique_plip_interactions"] = df["hash"].apply(lambda x: len(set(str(x).split(";"))) if x is not None and str(x) != "nan" else 0)
    df["has_pocket"] = df["single_pocket_ID"].apply(lambda x: x in cif_data["pockets"])
    df["UniProt_IDs"] = df.apply(lambda row: [cif_data["uniprot_ids"].get(f'{row["PDB_ID"]}_{y}', '') for y in sorted(row["prox_plip_chains"])] if str(row["prox_plip_chains"]) != "nan" else np.nan, axis=1)
    uniprot_to_chains = get_uniprot_to_chains(cif_data["uniprot_ids"])
    df["num_pdbs_for_uniprots"] = df["UniProt_IDs"].apply(lambda x: [len(uniprot_to_chains.get(y, [])) for y in x] if str(x) != "nan" else np.nan)
    df["date"] = df["PDB_ID"].apply(lambda x: cif_data["dates"].get(x.upper(), np.nan))
    df["year"] = df["date"].apply(lambda x: x.split("-")[0] if str(x) != "nan" else np.nan)
    df["single_pocket_chains"] = df.apply(lambda row: "_".join(sorted(row["prox_plip_chains"])) if str(row["prox_plip_chains"]) != "nan" else np.nan, axis=1)
    df["protein_chain_lengths"] = df.apply(lambda row: [cif_data["lengths"][row["PDB_ID"]].get(y, np.nan) for y in row["single_pocket_chains"].split("_")] if row["PDB_ID"] in cif_data["lengths"] and str(row["single_pocket_chains"]) != "nan" else np.nan, axis=1)
    df["protein_entity_ids"] = df.apply(lambda row: [cif_data["entity_mapping"][row["PDB_ID"]].get(y, np.nan) for y in row["single_pocket_chains"].split("_")] if row["PDB_ID"] in cif_data["entity_mapping"] and str(row["single_pocket_chains"]) != "nan" else np.nan, axis=1)
    df["ligand_entity_id"] = df.apply(lambda row: cif_data["entity_mapping"][row["PDB_ID"]].get(row["ligand_mmcif_chain"], np.nan) if row["PDB_ID"] in cif_data["entity_mapping"] else np.nan, axis=1)
    return df


def assign_pocket(df):
    df["Pocket_Number"] = [None] * len(df)
    # Assign to same pocket if prox_ligand_chains overlaps with ligand_mmcif_chain
    for (pdb_id, biounit), group in tqdm(df.groupby(['PDB_ID', 'biounit'])):
        pocket_number = 0
        pocket_dict = {}
        for index, row in group.iterrows():
            prox_ligand_chains = row['prox_plip_ligand_chains']
            if prox_ligand_chains is None:
                prox_ligand_chains = set()
            prox_ligand_chains = set(prox_ligand_chains)
            assigned = False
            for key, value in pocket_dict.items():
                if row['ligand_mmcif_chain'] in value or len(prox_ligand_chains.intersection(value)) > 0:
                    df.loc[index, 'Pocket_Number'] = key
                    pocket_dict[key] |= prox_ligand_chains
                    assigned = True
                    break
            if not assigned:
                pocket_name = f"{pdb_id}_{int(biounit)}_{pocket_number}"
                ligand_chains = prox_ligand_chains.copy()
                ligand_chains.add(row['ligand_mmcif_chain'])
                pocket_dict[pocket_name] = ligand_chains
                df.loc[index, 'Pocket_Number'] = pocket_name
                pocket_number += 1
    df = df[df["Pocket_Number"].notnull()].reset_index(drop=True)
    pocket_to_id = {}
    pocket_to_nums = {}
    for _, pocket in tqdm(df.groupby("Pocket_Number")):
        ligand_chain_to_name = {}
        for _, row in pocket.iterrows():
            ligand_chain_to_name[row["ligand_mmcif_chain"]] = row["Ligand"]
        ligand_chains = sorted(ligand_chain_to_name)
        ligand_names = [ligand_chain_to_name[x] for x in sorted(ligand_chain_to_name)]
        protein_chains = sorted(set(y for x in pocket["prox_plip_chains"] for y in x if str(y) != "nan"))
        pocket_to_id[pocket["Pocket_Number"].values[0]] = f"{pocket['PDB_ID'].values[0]}__{int(pocket['biounit'].values[0])}__{'_'.join(protein_chains)}__{'_'.join(ligand_chains)}__{'_'.join(ligand_names)}"
        pocket_to_nums[pocket["Pocket_Number"].values[0]] = (len(ligand_chains), len(protein_chains))
    df["pocket_ID"] = df["Pocket_Number"].apply(lambda x: pocket_to_id[x])
    df["num_ligands_in_pocket"] = df["Pocket_Number"].apply(lambda x: pocket_to_nums[x][0])
    df["num_chains_in_pocket"] = df["Pocket_Number"].apply(lambda x: pocket_to_nums[x][1])
    return df

def get_smtl_info(pdb_id, smtl_dir):
    #TODO: get these from cif_data instead
    smtl_dir = Path(smtl_dir)
    pdb_id = pdb_id.lower()
    data = []
    for annotation_file in (smtl_dir / pdb_id[:2] / pdb_id[2:]).glob("annotation.*.json"):
        with open(annotation_file) as f:
            annotation = json.load(f)
        if annotation["status"] == "deleted":
            continue
        if not all(x.isdigit() for x in annotation["mmcif_id"]):
            continue
        data.append(dict(pdb_id=pdb_id.upper(), biounit=int(annotation["mmcif_id"]), 
                         method=annotation["method"],
                    resolution=annotation["resolution"], oligo_state=annotation["oligo_state"],
                    transmembrane=annotation.get("membrane", {}).get("is_transmem", None)))
    return data

def get_covalent_info(pdb_id, pdb_dir):
    from mmcif.io.PdbxReader import PdbxReader
    import gzip
    cif_file = Path(pdb_dir) / pdb_id[1:3].lower() / f"{pdb_id.lower()}.cif.gz"
    if not cif_file.exists():
        return None
    with gzip.open(str(cif_file), 'rt', encoding='utf-8') as f:
        data = []
        try:
            prd = PdbxReader(f)
            prd.read(data)
        except Exception as e:
            print(f"Error reading {cif_file}: {e}")
            return None
    struct_conn = data[0].getObj('struct_conn')
    if struct_conn is None:
        return None
    covalent_bonds = [index for index in range(struct_conn.getRowCount()) if struct_conn.getValue('conn_type_id', index) == 'covale']
    if len(covalent_bonds) == 0:
        return None
    covalent_ligand_chains = set()
    for index in covalent_bonds:
        covalent_ligand_chains.add(struct_conn.getValue('ptnr1_label_asym_id', index))
        covalent_ligand_chains.add(struct_conn.getValue('ptnr2_label_asym_id', index))
    return {pdb_id: covalent_ligand_chains}

def label_smtl(df, smtl_dir, num_threads=20):
    annotations = {}
    input_args = [(pdb_id, smtl_dir) for pdb_id in df["PDB_ID"].unique()]
    with Pool(num_threads) as p:
        for annotation_list in p.starmap(get_smtl_info, input_args):
            for annotation in annotation_list:
                annotations[(annotation["pdb_id"], annotation["biounit"])] = annotation
    for a in ["method", "oligo_state", "transmembrane", "resolution"]:
        df[a] = df.apply(lambda row: annotations.get((row["PDB_ID"], int(row["biounit"])), {}).get(a, np.nan), axis=1)
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
            for pdb_id, rows in tqdm(df.groupby("PDB_ID")):
                asym_auth_chain_mapping = {v: k for k, v in chain_mapping.get(pdb_id, {}).items()}
                pdb_id = pdb_id.lower()
                chains = set()
                for _, row in rows.iterrows():
                    chains |= set(row["prox_plip_chains"])
                auth_chains = set(asym_auth_chain_mapping[x] for x in chains if x in asym_auth_chain_mapping)
                for chain in auth_chains:
                    if source == "foldseek":
                        f.write(f"{pdb_id}.cif.gz_{chain}\n")
                    else:
                        f.write(f"{pdb_id}_{chain}\n")
                prefix = pdb_id[1:3]
                if prefix not in first_time:
                    first_time.add(prefix)
                    open_mode = "w"
                else:
                    open_mode = "a"
                with open(subset_folder / f"{prefix}.txt", open_mode) as fp:
                    for chain in auth_chains:
                        if source == "foldseek":
                            fp.write(f"{pdb_id}.cif.gz_{chain}\n")
                        else:
                            fp.write(f"{pdb_id}_{chain}\n")

def create_dataset_files(dataset_dir, foldseek_dir, mmseqs_dir, plip_file, validation_file, cif_data_dir, components_file, cofactor_file, artifact_file, 
                         smtl_dir, num_threads=20, max_protein_chains=10, max_ligand_chains=5, overwrite=False):
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
    smtl_dir: directory with SMTL files
    max_protein_chains: maximum number of protein chains to consider
    max_ligand_chains: maximum number of ligand chains to consider
    num_threads: number of threads to use
    overwrite: whether to overwrite existing files
    """
    dataset_dir = Path(dataset_dir)
    dataset_dir.mkdir(exist_ok=True)
    all_pockets_file = dataset_dir / "all_pockets.csv"
    
    if overwrite or not all_pockets_file.exists():
        plip_df = read_df(plip_file)
        plip_df = plip_df[(~plip_df["biounit"].str.startswith("ASU")) & (~plip_df["biounit"].str.startswith("PAU")) & (~plip_df["biounit"].str.startswith("XAU"))].reset_index(drop=True)
        print("Number of PLIP pockets:", len(plip_df))
        df = read_df(validation_file, columns=VALIDATION_COLUMNS)
        print("Number of validation pockets:", len(df))
        df.rename(columns={"PocketID": "validation_pocket_ID"}, inplace=True)
        df = pnd.merge(df, plip_df, 
                    on=['joint_pocket_ID', "PDB_ID", "Ligand", "ligand_mmcif_chain"], 
                    how='outer')
        df["biounit"] = df["biounit"].fillna(1)
        df["single_pocket_ID"] = df.apply(lambda row: f'{row["PDB_ID"]}__{int(row["biounit"])}__{row["ligand_mmcif_chain"]}', axis=1)
        print("Annotation of chains and residues")
        df = annotate_chains_and_residues(df, cif_data_dir)
        RDLogger.DisableLog('rdApp.*')
        parsed_components = ccd_reader.read_pdb_components_file(str(components_file))
        print("Labeling ligand types")
        df = label_ligand_types(df, parsed_components, cofactor_file)
        print("Labeling artifacts")
        df = label_artifacts(df, artifact_file)
        print("Annotating info from SMTL")
        df = label_smtl(df, smtl_dir, num_threads=num_threads)
        df["has_plip"] = df["plip_pocket_ID"].notna()
        df["has_validation"] = df["validation_pocket_ID"].notna()
        df.to_csv(all_pockets_file, index=False)
    df = read_df(all_pockets_file)
    print("Number of pockets after merging PLIP and validation:", len(df), df["joint_pocket_ID"].nunique())
    df = df[df["has_plip"] & df["has_pocket"]].reset_index(drop=True)
    df = assign_pocket(df)
    print("Number of pockets after filtering for PLIP and binding sites and merging pockets")
    print("\tTotal rows:", len(df))
    print("\tTotal pockets:", df["joint_pocket_ID"].nunique())
    print("\tTotal merged pockets:", df["Pocket_Number"].nunique())
    df.to_csv(dataset_dir / "all_pockets_with_plip_and_bs.csv", index=False)
    df["is_ion"] = df["ligand_type"] == "ion"
    only_ions = df.groupby("Pocket_Number").agg({"is_ion": "all"})
    only_ions = only_ions[only_ions["is_ion"]].index
    df = df[~df["Pocket_Number"].isin(only_ions)]
    df.to_csv(dataset_dir / "small_molecule_pockets.csv", index=False)
    print("Number of pockets with small molecules:", df["Pocket_Number"].nunique())
    sm_pockets = df[(~df["is_ion"]) & (~df["is_artifact"])]["Pocket_Number"].unique()
    df = df[df["Pocket_Number"].isin(sm_pockets)]
    print("Number of pockets with small molecules (no artifacts):", df["Pocket_Number"].nunique())
    df.to_csv(dataset_dir / "small_molecule_pockets_no_artifacts.csv", index=False)
    df = df[(df["num_chains_in_pocket"] <= max_protein_chains) & (df["num_ligands_in_pocket"] <= max_ligand_chains)]
    print("Number of pockets with small molecules (no artifacts, max protein chains, max ligand chains):", df["Pocket_Number"].nunique())
    df.to_csv(dataset_dir / "filtered_pockets.csv", index=False)
    chain_mapping = load_cif_data(cif_data_dir, names_to_load=["chain_mapping"])["chain_mapping"]
    make_subset_files(df, Path(foldseek_dir) / "filtered_pockets", Path(mmseqs_dir) / "filtered_pockets", chain_mapping, "filtered_pockets")

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
    parser.add_argument("smtl_dir", type=str, help="Directory with SMTL files")
    parser.add_argument("--num_threads", type=int, default=20, help="Number of threads to use")
    parser.add_argument("--max_protein_chains", type=int, default=10, help="Maximum number of protein chains to consider")
    parser.add_argument("--max_ligand_chains", type=int, default=5, help="Maximum number of ligand chains to consider")
    parser.add_argument("--overwrite", action="store_true", help="Whether to overwrite existing files")

    args = parser.parse_args()
    create_dataset_files(**vars(args))

if __name__ == "__main__":
    main()