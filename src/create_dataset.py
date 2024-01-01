from multiprocessing import Pool
import pandas as pnd
from tqdm import tqdm
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pnd
from collections import defaultdict
from matplotlib.colors import hex2color
import json
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*') 
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
    for column in ["center_of_mass", "protein_chain_lengths", "pdb_id_protein_chains", "UniProt_IDs", "num_uniprot_ids", "num_pdbs_for_uniprots"]:
        if column in df.columns:
            df[column] = df[column].apply(lambda x: parse_list_set(x))
    return df


def parse_cofactors(cofactor_file):
    """Cofactors from PDBe and https://doi.org/10.1021/acs.jcim.3c01573"""
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
    df["is_artifact"] = df["is_artifact"] | (df["num_prox_plip_residues"] < num_prox_plip_residues_threshold)
    ligand_counts = dict(zip(df["Ligand"].value_counts().index, df["Ligand"].value_counts().values))
    common_ligands = set(str(x) for x in ligand_counts if ligand_counts[x] > common_ligand_threshold)
    df["is_artifact"] = df["is_artifact"] | ((df["Ligand"].isin(ignore_ligands)) & (df["Ligand"].isin(common_ligands)) & (df["ligand_type"] != "cofactor"))
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


def annotate_chains_and_residues(df, pocket_dir):
    pocket_residues = dict()
    for bs_file in tqdm(Path(pocket_dir).iterdir()):
        with open(bs_file) as f:
            pocket_residues.update(json.load(f))
    df["num_prox_plip_chains"] = df["prox_plip_chains"].apply(lambda x: len(set(x)) if str(x) != "nan" else 0)
    df["num_prox_plip_residues"] = df["prox_plip_residues"].apply(lambda x: len(set(x)) if str(x) != "nan" else 0)
    df["prox_chains"] = df["joint_pocket_ID"].apply(lambda x: set(y["chain"] for y in pocket_residues[x]["binding_site"]) if x in pocket_residues else set())
    df["num_prox_chains"] = df["prox_chains"].apply(lambda x: len(x))
    df["prox_residues"] = df["joint_pocket_ID"].apply(lambda x: set(f'{y["chain"]}:{y["residue_number"]}' for y in pocket_residues[x]["binding_site"]) if x in pocket_residues else set())
    df["center_of_mass"] = df["joint_pocket_ID"].apply(lambda x: pocket_residues[x]["center_of_mass"] if x in pocket_residues else None)
    df["num_plip_interactions"] = df["hash"].apply(lambda x: len(str(x).split(";")) if x is not None and str(x) != "nan" else 0)
    df["num_unique_plip_interactions"] = df["hash"].apply(lambda x: len(set(str(x).split(";"))) if x is not None and str(x) != "nan" else 0)
    df["has_bs"] = df["joint_pocket_ID"].apply(lambda x: x in pocket_residues)
    return df


def assign_pocket(df):
    df["Pocket_Number"] = [None] * len(df)
    # Assign to same pocket if prox_ligand_chains overlaps with ligand_mmcif_chain
    for (pdb_id, biounit), group in tqdm(df.groupby(['PDB_ID', 'biounit'])):
        pocket_number = 0
        pocket_dict = {}
        for index, row in group.iterrows():
            assigned = False
            for key, value in pocket_dict.items():
                if row['ligand_mmcif_chain'] in value or len(set(row['prox_plip_ligand_chains']).intersection(value)) > 0:
                    df.loc[index, 'Pocket_Number'] = key
                    pocket_dict[key] |= row['prox_plip_ligand_chains']
                    assigned = True
                    break
            if not assigned:
                pocket_name = f"{pdb_id}_{int(biounit)}_{pocket_number}"
                ligand_chains = row['prox_plip_ligand_chains']
                ligand_chains.add(row['ligand_mmcif_chain'])
                pocket_dict[pocket_name] = ligand_chains
                df.loc[index, 'Pocket_Number'] = pocket_name
                pocket_number += 1
    df = df[df["Pocket_Number"].notnull()].reset_index(drop=True)
    number_to_name = df.groupby("Pocket_Number").apply(lambda x: "_".join(sorted(set(x["ligand_mmcif_chain"])))).to_dict()
    df["Pocket_Name"] = df["Pocket_Number"].apply(lambda x: f"{x.split('_')[0]}__{x.split('_')[1]}__{number_to_name[x]}" if str(x) != "nan" else np.nan)
    df["num_ligands_in_pocket"] = df["Pocket_Name"].apply(lambda x: len(x.split("__")[-1].split("_")) if str(x) != "nan" else 0)
    return df
    

def label_smtl_info(pdb_id, biounit, smtl_dir):
    pdb_id = pdb_id.lower()
    annotation_file = smtl_dir / pdb_id[:2] / pdb_id[2:] / f"annotation.{biounit}.json"
    if not annotation_file.exists():
        return {}
    with open(annotation_file) as f:
        annotation = json.load(f)
    lengths = {}
    for x in annotation['entities']:
        if "seqres" in x:
            for c in x['chains']:
                if c['orig_pdb_name'] is not None:
                    key = f"{pdb_id}_{c['orig_cif_name']}"
                    lengths[key] = len(x["seqres"])
    annotation = {"lengths": lengths, "method": annotation["method"], 
                    "resolution": annotation["resolution"],
                  "oligo_state": annotation["oligo_state"], 
                  "transmembrane": annotation.get("membrane", {}).get("is_transmem", None)}
    return annotation

def label_uniprots(df, pdb_chain_uniprot_mapping_file, uniprot_pdb_mapping_file):
    pdb_uniprot = pnd.read_csv(pdb_chain_uniprot_mapping_file, sep="\t", comment="#")
    pdb_uniprot["pdb_id_protein_chain"] = pdb_uniprot["PDB"] + "_" + pdb_uniprot["CHAIN"]
    pdb_chain_to_uniprots = defaultdict(set)
    for _, row in tqdm(pdb_uniprot.iterrows()):
        pdb_chain_to_uniprots[row["pdb_id_protein_chain"]].add(row["SP_PRIMARY"])
    df["UniProt_IDs"] = df["pdb_id_protein_chains"].apply(lambda x: [";".join(pdb_chain_to_uniprots.get(y, set())) for y in x])
    df["num_uniprot_ids"] = df["UniProt_IDs"].apply(lambda x: [len(y.split(";")) for y in x])
    uniprot_df = pnd.read_csv(uniprot_pdb_mapping_file, comment="#", sep="\t")
    uniprot_df["num_pdbs"] = uniprot_df["PDB"].apply(lambda x: len(x.split(";")))
    id_to_num_pdbs = dict(zip(uniprot_df["SP_PRIMARY"], uniprot_df["num_pdbs"]))
    df["num_pdbs_for_uniprots"] = df["UniProt_IDs"].apply(lambda x: [sum(id_to_num_pdbs.get(z, 0) for z in y.split(";")) for y in x] if str(x) != "nan" else 0)
    return df

def label_dates(df, date_dir):
    dates = dict()
    for date_file in tqdm(Path(date_dir).iterdir()):
        with open(date_file) as f:
            dates.update(json.load(f))
    df["date"] = df["PDB_ID"].apply(lambda x: dates.get(x.upper(), np.nan))
    df["year"] = df["date"].apply(lambda x: x.split("-")[0] if str(x) != "nan" else np.nan)
    return df

def annotate_df(df, smtl_dir, date_dir, pdb_chain_uniprot_mapping_file, uniprot_pdb_mapping_file, num_threads=20):
    df["pdb_id_protein_chains"] = df.apply(lambda x: [f"{x['PDB_ID'].lower()}_{y}" for y in sorted(x['prox_plip_chains'])] if str(x['prox_plip_chains']) != "nan" else [], axis=1)
    annotations = {}
    input_args = [(row["PDB_ID"], row["biounit"], smtl_dir) for _, row in df[["PDB_ID", "biounit"]].drop_duplicates().iterrows()]
    with Pool(num_threads) as p:
        for (pdb_id, biounit, _), annotation in zip(p.starmap(label_smtl_info, input_args), input_args):
            annotations[(pdb_id, biounit)] = annotation
    for a in ["method", "oligo_state", "transmembrane", "resolution"]:
        df[a] = df.apply(lambda row: annotations[(row["PDB_ID"], row["biounit"])].get(a, np.nan), axis=1)
    df["protein_chain_lengths"] = df["pdb_id_protein_chains"].apply(lambda x: [annotations[(y.split("_")[0].upper(), 1)]["lengths"].get(y, np.nan) for y in x])
    df = label_dates(df, date_dir)
    df = label_uniprots(df, pdb_chain_uniprot_mapping_file, uniprot_pdb_mapping_file)
    df["single_pocket_ID"] = df.apply(lambda row: f'{row["PDB_ID"]}_{int(row["biounit"])}__{"_".join(sorted(row["prox_plip_chains"]))}__{row["ligand_mmcif_chain"]}__{row["Ligand"]}' if str(row["prox_plip_chains"]) != "nan" else np.nan, axis=1)
    return df

def create_dataset_files(dataset_dir, plip_file, validation_file, pocket_dir, components_file, cofactor_file, artifact_file, 
                         smtl_dir, date_dir, pdb_chain_uniprot_mapping_file, uniprot_pdb_mapping_file, num_threads=20, overwrite=False):
    """
    dataset_dir: directory to save the dataset files
    plip_file: file with PLIP results
    validation_file: file with validation results
    pocket_dir: directory with binding site files
    components_file: file with PDB chemical component dictionary
    cofactor_file: file with cofactors (https://www.ebi.ac.uk/pdbe/api/pdb/compound/cofactors)
    artifact_file: file with artifacts (https://github.com/kad-ecoli/mmCIF2BioLiP/blob/dc9769f286eafc550f799239486ef64450728246/ligand_list)
    smtl_dir: directory with SMTL files
    date_dir: directory with PDB dates
    pdb_chain_uniprot_mapping_file: file with PDB to UniProt mapping (ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz)
    uniprot_pdb_mapping_file: file with UniProt to PDB mapping (ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/uniprot_pdb.tsv.gz)
    num_threads: number of threads to use
    overwrite: whether to overwrite existing files
    """
    dataset_dir = Path(dataset_dir)
    dataset_dir.mkdir(exist_ok=True)
    all_pockets_file = dataset_dir / "all_pockets.csv"
    if overwrite or not all_pockets_file.exists():
        plip_df = read_df(plip_file)
        print("Number of PLIP pockets:", len(plip_df))
        df = read_df(validation_file, columns=VALIDATION_COLUMNS)
        print("Number of validation pockets:", len(df))
        df.rename(columns={"PocketID": "validation_pocket_ID"}, inplace=True)
        df = pnd.merge(df, plip_df, 
                    on=['joint_pocket_ID', "PDB_ID", "Ligand", "ligand_mmcif_chain"], 
                    how='outer')
        df = annotate_chains_and_residues(df, pocket_dir)
        RDLogger.DisableLog('rdApp.*')
        parsed_components = ccd_reader.read_pdb_components_file(str(components_file))
        df = label_ligand_types(df, parsed_components, cofactor_file)
        df = label_artifacts(df, artifact_file)
        df = annotate_df(df, smtl_dir, date_dir, pdb_chain_uniprot_mapping_file, uniprot_pdb_mapping_file, num_threads=num_threads)
        df["has_plip"] = df["plip_pocket_ID"].notna()
        df["has_validation"] = df["validation_pocket_ID"].notna()
        df.to_csv(all_pockets_file, index=False)
    df = read_df(all_pockets_file)
    print("Number of pockets after merging PLIP and validation:", len(df), df["joint_pocket_ID"].nunique())
    df = assign_pocket(df)
    print("Number of pockets after assigning pockets:", len(df), df["joint_pocket_ID"].nunique(), df["Pocket_Number"].nunique())
    df = df[df["has_plip"] & df["has_bs"]]
    print("Number of pockets after filtering for PLIP and binding sites")
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

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset_dir", type=str, help="Directory to save the dataset files")
    parser.add_argument("plip_file", type=str, help="File with PLIP results")
    parser.add_argument("validation_file", type=str, help="File with validation results")
    parser.add_argument("pocket_dir", type=str, help="Directory with binding site files")
    parser.add_argument("components_file", type=str, help="File with PDB chemical component dictionary")
    parser.add_argument("cofactor_file", type=str, help="File with cofactors (https://www.ebi.ac.uk/pdbe/api/pdb/compound/cofactors)")
    parser.add_argument("artifact_file", type=str, help="File with artifacts (https://github.com/kad-ecoli/mmCIF2BioLiP/blob/dc9769f286eafc550f799239486ef64450728246/ligand_list)")
    parser.add_argument("smtl_dir", type=str, help="Directory with SMTL files")
    parser.add_argument("date_dir", type=str, help="Directory with PDB dates")
    parser.add_argument("pdb_chain_uniprot_mapping_file", type=str, help="File with PDB to UniProt mapping (ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz)")
    parser.add_argument("uniprot_pdb_mapping_file", type=str, help="File with UniProt to PDB mapping (ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/uniprot_pdb.tsv.gz)")
    parser.add_argument("--num_threads", type=int, default=20, help="Number of threads to use")

    args = parser.parse_args()
    create_dataset_files(**vars(args))

if __name__ == "__main__":
    main()