from dataclasses import dataclass, field
import dataclasses
import typing as ty
from pathlib import Path
import pandas as pnd
import json
from multiprocessing import Pool
from tqdm import tqdm
from collections import defaultdict

@dataclass 
class PLIPHash:
    chain: str
    residue: int
    type: str
    three_letter_code: str
    attributes: ty.Dict[str, str]

    @classmethod
    def from_row(cls, row):
        attributes = {}
        if row["type"] == "hydrogen_bonds":
            attributes = {"donortype": row["info_donortype"], "acceptortype": row["info_acceptortype"], "protisdon": row["info_protisdon"], "sidechain": row["info_sidechain"]}
        elif row["type"] == "hydrophobic_interactions":
            attributes = {}
        elif row["type"] == "salt_bridges":
            attributes = {"lig_group": row["info_lig_group"], "protispos": row["info_protispos"]}
        elif row["type"] == "water_bridges":
            attributes = {"donortype": row["info_donortype"], "acceptortype": row["info_acceptortype"], "protisdon": row["info_protisdon"]}
        elif row["type"] == "metal_complexes":
            attributes = {"metal_type": row["info_metal_type"], "target_type": row["info_target_type"], "coordination": row["info_coordination"], "geometry": row["info_geometry"], "location": row["info_location"]}
        elif row["type"] == "pi_stacks":
            attributes = {"type": row["info_type"]}
        elif row["type"] == "pi_cation_interactions":
            attributes = {"lig_group": row["info_lig_group"], "protcharged": row["info_protcharged"]}
        elif row["type"] == "halogen_bonds":
            attributes = {"donortype": row["info_donortype"], "acceptortype": row["info_acceptortype"]}
        return cls(row["orig_cif_chain"], row["position"], row["type"], row["three_letters"], attributes)
    
    def to_string(self, hash_type="all"):
        if hash_type == "all": return "__".join(str(x) for x in [self.chain, self.residue, self.type, self.three_letter_code]) + "__" + "__".join(f"{k}:{v}" for k, v in self.attributes.items())
        elif hash_type == "_residue": return "__".join(str(x) for x in [self.type, self.three_letter_code]) + "__" + "__".join(f"{k}:{v}" for k, v in self.attributes.items())
        elif hash_type == "_nowater": 
            if self.type == "water_bridges": 
                return None 
            else: 
                return self.to_string(hash_type="")
        elif hash_type == "": return f"{self.type}__" + "__".join(f"{k}:{v}" for k, v in self.attributes.items())
        else: raise ValueError(f"Invalid plip type: {hash_type}")

    @classmethod
    def from_string(cls, string):
        chain, residue, type, three_letters, *attributes = string.split("__")
        if ':' in attributes[0]:
            attributes = dict(x.split(":") for x in attributes)
        else:
            attributes = {}
        return cls(chain, residue, type, three_letters, attributes)

@dataclass
class PLIPInteraction:
    pdb_id: str
    biounit: str
    orig_pdb_chain: str
    orig_cif_chain: str
    orig_pdb_chain_ligand: str
    orig_cif_chain_ligand: str
    prox_ligand_chains: ty.Set[str]
    bond_id: str
    chain: str
    type: str
    htype: str
    distance: float
    position: int
    three_letters: str
    one_letter: str
    ligcoo: ty.List[float]
    protcoo: ty.List[float]
    ligtype: ty.Optional[str] = None
    extra: ty.Dict[str, ty.Any] = field(default_factory=dict)

    def __init__(self, pdb_id, biounit, orig_pdb_chain, orig_cif_chain, orig_pdb_chain_ligand, 
                 orig_cif_chain_ligand, prox_ligand_chains, bond_id, chain, type, htype, distance, position, three_letters, 
                 one_letter, ligcoo, protcoo, ligtype, **kwargs):
        self.pdb_id = pdb_id
        self.biounit = biounit
        self.orig_pdb_chain = orig_pdb_chain
        self.orig_cif_chain = orig_cif_chain
        self.orig_pdb_chain_ligand = orig_pdb_chain_ligand
        self.orig_cif_chain_ligand = orig_cif_chain_ligand
        self.prox_ligand_chains = prox_ligand_chains
        self.bond_id = bond_id
        self.chain = chain
        self.type = type
        self.htype = htype
        self.distance = distance
        self.position = position
        self.three_letters = three_letters
        self.one_letter = one_letter
        self.ligcoo = ligcoo
        self.protcoo = protcoo
        self.ligtype = ligtype
        self.extra = kwargs

    @classmethod
    def from_json(cls, pdb_id, biounit, orig_pdb_chain, orig_cif_chain, orig_pdb_chain_ligand, orig_cif_chain_ligand, prox_ligand_chains, json_data):
        return cls(pdb_id, biounit, orig_pdb_chain, orig_cif_chain, orig_pdb_chain_ligand, orig_cif_chain_ligand, prox_ligand_chains, **json_data)


def interactions_to_dataframe(interactions: ty.List[PLIPInteraction]) -> pnd.DataFrame:
    """Converts a list of PLIPInteraction instances into a pandas DataFrame."""
    df = pnd.DataFrame([dataclasses.asdict(interaction) for interaction in interactions])
    if df.empty:
        return df
    df = df.join(df["extra"].apply(pnd.Series).add_prefix("info_"))
    df = df.join(df["info_extra"].apply(pnd.Series).add_prefix("info_"))
    df = df.drop(columns=["extra", "info_extra"])
    df["hash"] = df.apply(lambda row: PLIPHash.from_row(row).to_string(), axis=1)
    df["info_restype_lig"] = df["info_restype_lig"].apply(lambda x: "SODIUM" if str(x) == "nan" else x)
    df = df[df["hash"].notnull()].reset_index(drop=True)
    df = df.drop_duplicates(subset=["pdb_id", "biounit", "orig_cif_chain", "orig_cif_chain_ligand", "position", "hash"], keep="first")
    return df

def get_interactions_from_json(pdb_id, json_file: Path):
    """Extracts interactions from a PLIP JSON file."""
    with open(json_file) as f:
        json_data = json.load(f)
    if json_data["status"] == "deleted" or 'plip' not in json_data or json_data['plip'] is None:
        return [], []
    biounit = json_data["mmcif_id"]
    # if biounit is not an integer, return (skip asymmetric units)
    if not all(x.isdigit() for x in biounit):
        return [], []
    # Get the mapping between the original MMCIF chains and the chains used by PLIP
    orig_pdb_mapping = {}
    orig_cif_mapping = {}
    interacting_ligands = defaultdict(set)
    protein_chains = set()
    for x in json_data['entities']:
        for c in x['chains']:
            if c['orig_pdb_name'] is not None:
                if 'atom_seq' in c:
                    protein_chains.add(c['orig_cif_name'])
                orig_pdb_mapping[c['name']] = c['orig_pdb_name']
                orig_cif_mapping[c['name']] = c['orig_cif_name']

    # Get the interacting ligands for each ligand chain
    for x in json_data['entities']:
        for c in x['chains']:
            if c['orig_pdb_name'] is not None:
                for ligand in c.get('in_contact', []):
                    chain, name = ligand.split(".")
                    if chain == "_":
                        interacting_ligands[c['orig_cif_name']].add(orig_cif_mapping[name])
    interactions = []
    errors = []
    for key in json_data['plip']:
        for interaction in json_data['plip'][key]:
            if interaction['id'].split(":")[1] != "_":
                continue
            try:
                orig_pdb_chain_ligand = orig_pdb_mapping[str(interaction['position'])]
                orig_cif_chain_ligand = orig_cif_mapping[str(interaction['position'])]
                ligtype = interaction['ligtype']
                for chain in interaction['chains']:
                    if chain not in orig_pdb_mapping:
                        continue
                    orig_pdb_chain = orig_pdb_mapping[chain]
                    orig_cif_chain = orig_cif_mapping[chain]
                    if orig_cif_chain in protein_chains:
                        for bond_type in interaction['chains'][chain]['interactions']:
                            for x in interaction['chains'][chain]['interactions'][bond_type]:
                                x['ligtype'] = ligtype
                                interactions.append(PLIPInteraction.from_json(pdb_id, biounit, orig_pdb_chain, orig_cif_chain, 
                                                                            orig_pdb_chain_ligand, orig_cif_chain_ligand, 
                                                                            interacting_ligands[orig_cif_chain_ligand], x))
            except Exception as e:
                if e == KeyboardInterrupt:
                    raise e
                errors.append((pdb_id, biounit, interaction['position'], e))
    return interactions, errors
    

def interactions_from_folder(super_pdb_folder_output_folder, overwrite=True):
    super_pdb_folder, output_folder = super_pdb_folder_output_folder
    output_folder = Path(output_folder)
    super_pdb_folder = Path(super_pdb_folder)
    output_file = output_folder / f"{super_pdb_folder.name}.tsv"
    if not overwrite and output_file.exists():
        return pnd.read_csv(output_file, sep="\t", low_memory=False)
    interactions = []
    errors = []
    for pdb_folder in super_pdb_folder.iterdir():
        if pdb_folder.name.startswith(".") or not pdb_folder.is_dir():
            continue
        pdb_id = f"{super_pdb_folder.name}{pdb_folder.name}"
        for json_file in pdb_folder.iterdir():
            if not json_file.stem.startswith("annotation."):
                continue
            interactions_, errors_ = get_interactions_from_json(pdb_id, json_file)
            interactions += interactions_
            errors += errors_
    df = interactions_to_dataframe(interactions)
    df.to_csv(output_file, sep="\t", index=False)
    errors_file = output_folder / f"{super_pdb_folder.name}_errors.tsv"
    errors_file.unlink(missing_ok=True)
    if len(errors):
        with open(errors_file, "w") as f:
            f.write("pdb_id\tbiounit\tposition\terror\n")
            for error in errors:
                f.write("\t".join(str(x) for x in error) + "\n")
    return df

def group_interactions(output_folder):
    """
    Groups all interactions belonging to a single ligand pocket SPLC
    defined by the PDB ID, biounit, ligand name and ligand MMCIF chain.
    """
    output_folder = Path(output_folder)
    df_interactions = []
    for plip_file in tqdm(output_folder.iterdir()):
        if not plip_file.name.endswith(".tsv") or "errors" in plip_file.name or "interactions" in plip_file.name:
            continue
        try:
            df = pnd.read_csv(plip_file, sep="\t", low_memory=False)
            df["prox_ligand_chains"] = df["prox_ligand_chains"].apply(lambda x: eval(x))
        except pnd.errors.EmptyDataError:
            continue
        if df.empty:
            continue
        df["orig_cif_chain"] = df["orig_cif_chain"].apply(lambda x: "NA" if str(x) == "nan" else x)
        df["orig_cif_chain_ligand"] = df["orig_cif_chain_ligand"].apply(lambda x: "NA" if str(x) == "nan" else x)
        df["info_restype_lig"] = df["info_restype_lig"].apply(lambda x: "SODIUM" if str(x) == "nan" else x)
        df["plip_pocket_ID"] = ["_".join(x) for x in zip(df["pdb_id"].str.upper(), df["biounit"].astype(str), df["info_restype_lig"].astype(str), df["orig_cif_chain_ligand"].astype(str))]       
        df["protein_residue"] = df["orig_cif_chain"] + ":" + df["position"].astype(str)
        df_prox_res = df.groupby("plip_pocket_ID").agg({"protein_residue": "unique"})
        df_prox_chains = df.groupby("plip_pocket_ID").agg({"orig_cif_chain": "unique"})
        df_prox_ligand_chains = df.groupby("plip_pocket_ID").agg({"prox_ligand_chains": "first"})
        plip_df = df.groupby(["plip_pocket_ID", "ligtype"])["hash"].apply(lambda x: ";".join(sorted(x))).reset_index()
        plip_df["PDB_ID"] = plip_df["plip_pocket_ID"].apply(lambda x: x.split("_")[0].upper())
        plip_df["biounit"] = plip_df["plip_pocket_ID"].apply(lambda x: x.split("_")[1])
        plip_df["Ligand"] = plip_df["plip_pocket_ID"].apply(lambda x: x.split("_")[2].upper())
        plip_df["ligand_mmcif_chain"] = plip_df["plip_pocket_ID"].apply(lambda x: x.split("_")[3])
        plip_df["ligand_mmcif_chain"] = plip_df["ligand_mmcif_chain"].apply(lambda x: "NA" if str(x) == "nan" else x)
        plip_df["prox_plip_residues"] = plip_df["plip_pocket_ID"].map(df_prox_res["protein_residue"])
        plip_df["prox_plip_chains"] = plip_df["plip_pocket_ID"].map(df_prox_chains["orig_cif_chain"])
        plip_df["prox_plip_ligand_chains"] = plip_df["plip_pocket_ID"].map(df_prox_ligand_chains["prox_ligand_chains"])
        df_interactions.append(plip_df)
    df_interactions = pnd.concat(df_interactions, ignore_index=True)
    df_interactions.to_csv(output_folder / "interactions.tsv", sep="\t", index=False)


def main():
    from sys import argv
    smtl_dir, output_folder, n_threads = argv[1:]
    output_folder = Path(output_folder)
    smtl_dir = Path(smtl_dir)
    num_per_folder = {d: sum(1 for _ in d.iterdir()) for d in smtl_dir.iterdir()}
    sorted_names = sorted(num_per_folder, key=lambda x: num_per_folder[x], reverse=True)
    with Pool(int(n_threads)) as p:
        _ = list(tqdm(p.imap(interactions_from_folder, 
                             [(s, output_folder) for s in sorted_names]),
                             total=len(sorted_names)))
    group_interactions(output_folder)


if __name__ == "__main__":
    main()
    