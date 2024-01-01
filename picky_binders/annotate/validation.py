import numpy as np
import pandas as pnd
from tqdm import tqdm
from pathlib import Path

def is_not_float(x):
    try:
        x = float(x)
        return False
    except ValueError as e:
        return True
    
def is_nan(x):
    return x == "nan" or np.isnan(float(x))

THRESHOLDS = {
    "entry_resolution": 3.5,
    "entry_rfree": 0.45,
    "entry_r": 0.4,
    "entry_r_minus_rfree": 0.05,
    "rscc": 0.8,
    "rsr": 0.3,
    "avgoccu": 1.0,
    "perc_prox_rscc": 90
}

CRITERIA = {
    "entry_resolution": lambda x: float(x) <= THRESHOLDS["entry_resolution"],
    "entry_r": lambda x : is_not_float(x) or is_nan(x) or float(x) < THRESHOLDS["entry_r"],
    "entry_r_minus_rfree": lambda x : is_nan(x) or x <= THRESHOLDS["entry_r_minus_rfree"],
    "rscc": lambda x : float(x) > THRESHOLDS["rscc"],
    "rsr": lambda x : float(x) <= THRESHOLDS["rsr"],
    "altcode": lambda x: x == " ",
    "avgoccu": lambda x: float(x) == 1.0,
    "prox_alternative_configuration_residues_flag": lambda x: x == "False" or x == False,
    "perc_prox_rscc": lambda x: float(x) > THRESHOLDS["perc_prox_rscc"],
}

IRIDIUM_CRITERIA = {
    "entry_resolution": lambda x: float(x) <= 3.5,
    "entry_r": lambda x : is_not_float(x) or is_nan(x) or float(x) < 0.4,
    "entry_rfree": lambda x : is_not_float(x) or is_nan(x) or float(x) < 0.45,
    "entry_r_minus_rfree": lambda x : is_nan(x) or x <= 0.05,
    "rscc": lambda x : float(x) >= 0.9,
    "rsr": lambda x : float(x) <= 0.1,
    "avgoccu": lambda x: float(x) == 1.0,
    "altcode": lambda x: x == " ",
    "prox_alternative_configuration_residues_flag": lambda x: x == "False" or x == False,
    "prox_rscc": lambda x : all(float(y) >= 0.9 for y in str(x).split(";") if y != "nan"),
    "prox_rsr": lambda x : all(float(y) <= 0.1 for y in str(x).split(";") if y != "nan"),
    "prox_occupancies": lambda x: all(float(o) == 1.0 for y in str(x).replace("[","").replace("]","").split(";") for o in str(y).split(", ") ),
}

def pd_filter_iridium(row):
    for column in IRIDIUM_CRITERIA:
        if not IRIDIUM_CRITERIA[column](row[column]):
            return False
    return True

def pd_filter(row):
    for column in CRITERIA:
        if not CRITERIA[column](row[column]):
            return False
    return True

def filter_df(folder, parsed_components=None):
    """
    Filter the validation results by the criteria defined in the CRITERIA dictionary.
    """
    bad_files = []
    validation_df = []
    num_total_pockets = 0
    for filename in tqdm(folder.glob("Result*.csv")):
        df = pnd.read_csv(filename, dtype="object")
        num_total_pockets += len(df)
        if "entry_resolution" not in df.columns or "Unnamed: 0" not in df.columns:
            bad_files.append(filename.name)
            continue
        df = df[df["entry_resolution"].apply(lambda x: not str(x).startswith("molprobity"))]
        if len(df):
            df = df.rename(columns={"Unnamed: 0": "PocketID"})
            df["Ligand"] = df["PocketID"].map(lambda x: str(x).split(":")[-1].split(".")[0])
            df["PDB_chain"] = df["PocketID"].map(lambda x: str(x).split(":")[0])
            df["PDB_ID"] = df["PDB_chain"].map(lambda x: x[:4])
            validation_df.append(df)

    validation_df = pnd.concat(validation_df, ignore_index=True)
    print(f"Number of bad files: {len(bad_files)}")
    print(f"Number of total pockets: {num_total_pockets}")
    print(f"Number of bad files: {len(bad_files)}")
    print(f"Number of total pockets: {num_total_pockets}")
    validation_df = validation_df[~validation_df["NatomsEDS"].isna()]
    print(f"After filtering by EDS availability: {len(validation_df)} pockets across {len(validation_df['PDB_ID'].unique())} proteins.")

    # Mapping information
    validation_df["entry_r_minus_rfree"] = validation_df[["entry_r", "entry_rfree"]].apply(lambda x: np.abs(float(x[0]) - float(x[1])) if not is_not_float(x[0]) and not is_not_float(x[1]) else np.nan, axis=1)
    if parsed_components is not None:
        validation_df["ccd_NumRotatableBonds"] = validation_df["Ligand"].map(lambda x: parsed_components[x].component.physchem_properties.get("NumRotatableBonds", np.nan))
    validation_df["perc_prox_rscc"] = validation_df["prox_rscc"].map(lambda x: 100*sum(float(y) > THRESHOLDS["rscc"] for y in str(x).split(";")) / len(str(x).split(";")) if "nan" not in str(x) else np.nan)
    validation_df["iridium_pass"] = validation_df.apply(pd_filter_iridium, axis=1)
    validation_df["pass"] = validation_df.apply(pd_filter, axis=1)

def write_jobs(validation_super_folder, results_folder, logs_folder, script_path, commands_file):
    with open(commands_file, "w") as f:
        for folder in Path(validation_super_folder).iterdir():
            folder = folder.name
            f.write(f"{script_path} -v -n 1 {results_folder}/Results_{folder}.csv {logs_folder} {validation_super_folder}/{folder}/*/*_validation.xml.gz\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('validation_folder', type=str, help='Folder containing the validation results.')
    parser.add_argument('results_folder', type=str, help='Folder where the results will be saved.')
    parser.add_argument('logs_folder', type=str, help='Folder where the logs will be saved.')
    parser.add_argument('script_path', type=str, help='Path to the PDBValidationExtract script.')
    parser.add_argument('commands_file', type=str, help='Path to the file where the commands will be saved.')

    args = parser.parse_args()
    write_jobs(**vars(args))
