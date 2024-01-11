from ost import conop, mol, io
from pathlib import Path
from tqdm import tqdm
import pandas as pnd
from Bio.PDB import PDBParser
from Bio.PDB.mmcifio import MMCIFIO
from multiprocessing import Pool
import logging

# Define available names for protein and ligand chains
PROTEIN_CHAINS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
LIGAND_CHAINS = PROTEIN_CHAINS.lower() + '0123456789'

def filter_and_rename_chains(ent, protein_chains, ligand_chains):
    # Create a subset entity with the selected chains
    ent_selected = mol.CreateEntityFromView(ent.Select(" or ".join(f"chain='{c}'" for c in protein_chains + ligand_chains)), True)

    # Intermediate renaming step
    intermediate_names = {}
    edi = ent_selected.EditXCS(mol.BUFFERED_EDIT)
    for i, chain in enumerate(ent_selected.GetChainList()):
        intermediate_names[f"T{i}"] = chain.name
        edi.RenameChain(chain, f"T{i}")
    edi.UpdateICS()

    # Final renaming step
    protein_chain_index = 0
    ligand_chain_index = 0
    name_mapping = {}

    for chain in ent_selected.GetChainList():
        original_name = intermediate_names[chain.name]
        original_chain = ent.FindChain(original_name)
        if original_name in protein_chains:
            final_name = PROTEIN_CHAINS[protein_chain_index]
            protein_chain_index += 1
        elif original_name in ligand_chains:
            final_name = LIGAND_CHAINS[ligand_chain_index]
            ligand_chain_index += 1
        edi.RenameChain(chain, final_name)
        edi.SetChainDescription(chain, original_chain.description)
        edi.SetChainType(chain, original_chain.type)
        name_mapping[original_name] = final_name
    edi.UpdateICS()
    name_mapping_reverse = {v: k for k, v in name_mapping.items()}
    return name_mapping, name_mapping_reverse, ent_selected

def load_biounit_seqres(cif_file, biounit_id):
    try:
        ent, seqres, info = io.LoadMMCIF(str(cif_file), seqres=True, info=True, remote=False)
    except Exception as e:
        logging.error(f"Could not load {cif_file}: {e}")
        return None, None
    biounit_info = [b for b in info.biounits if b.id == biounit_id][0]
    try:
        biounit = mol.alg.CreateBU(ent, biounit_info)
    except Exception as e:
        logging.error(f"Could not create biounit for {cif_file}: {e}")
        return None, None
    return biounit, seqres

def save_ligands(new_ent, ligand_chains, output_prefix, name_mapping):
    for original_ligand_chain in ligand_chains:
        if original_ligand_chain not in name_mapping:
            logging.error(f"Could not find ligand chain {original_ligand_chain} in {output_prefix}")
            return
        new_ligand_chain = name_mapping[original_ligand_chain]
        ligand = new_ent.Select(f"chain='{new_ligand_chain}'")
        try:
            io.SaveEntity(ligand, f"{output_prefix}__{new_ligand_chain}__{original_ligand_chain}.sdf")
        except Exception as e:
            logging.error(f"Could not save ligand {original_ligand_chain} in {output_prefix}: {e}")
            return

def save_sequences(new_ent, output_fasta_file, chain_to_sequence, name_mapping_reverse):
     with open(output_fasta_file, "w") as f:
        for x in new_ent.chains:
            if x.name not in name_mapping_reverse:
                logging.error(f"Could not find chain {x.name} in {output_fasta_file}")
                return
            original_chain = name_mapping_reverse[x.name]
            if original_chain in chain_to_sequence:
                f.write(f">{x.name}__{original_chain}\n")
                f.write(chain_to_sequence[original_chain] + "\n")

def save_pdb_cif_files(new_ent, pocket_name, output_pdb_file, output_cif_file):
    #TODO replace biopython mmcif writer with ost MMCifWriter
    try:
        io.SavePDB(new_ent, str(output_pdb_file))
    except Exception as e:
        logging.error(f"Could not save {output_pdb_file}: {e}")
        return
    try:
        parser = PDBParser()
        structure = parser.get_structure(pocket_name, str(output_pdb_file))
        mmcif_io = MMCIFIO()
        mmcif_io.set_structure(structure)
        mmcif_io.save(str(output_cif_file))
    except Exception as e:
        logging.error(f"Could not save {output_cif_file}: {e}")
        return

def save_pocket(pocket_name, cif_folder, output_folder, compound_lib):
    pdb_id, biounit_id, protein_chains, ligand_chains, _ = pocket_name.split("__")
    protein_chains, ligand_chains = protein_chains.split("_"), ligand_chains.split("_")
    cif_file = Path(cif_folder) / f"{pdb_id.lower()[1:3]}" / f"{pdb_id.lower()}.cif.gz"
    if len(protein_chains) > len(PROTEIN_CHAINS) or len(ligand_chains) > len(LIGAND_CHAINS):
        logging.info(f"Too many chains for {pocket_name}")
        return
    if not cif_file.exists():
        logging.info(f"{cif_file} does not exist")
        return
    output_pdb_file = Path(output_folder) / "pdb_files" / f"{pocket_name}.pdb"
    output_cif_file = Path(output_folder) / "cif_files" / f"{pocket_name}.cif"
    output_fasta_file = Path(output_folder) / "fasta_files" / f"{pocket_name}.fasta"
    output_ligand_prefix = str(Path(output_folder) / "ligand_files" / f"{pocket_name}")
    clib = conop.CompoundLib.Load(compound_lib)
    conop.SetDefaultLib(clib) 
    
    biounit, seqres = load_biounit_seqres(cif_file, biounit_id)
    if biounit is None or seqres is None:
        return
    chain_to_sequence = {f"1.{s.name}": s.string for s in seqres if s.name in protein_chains}
    protein_chains = [f"1.{x}" for x in protein_chains]
    ligand_chains = [f"1.{x}" for x in ligand_chains]
    name_mapping, name_mapping_reverse, new_ent = filter_and_rename_chains(biounit, protein_chains, ligand_chains)
    save_ligands(new_ent, ligand_chains, output_ligand_prefix, name_mapping)
    save_sequences(new_ent, output_fasta_file, chain_to_sequence, name_mapping_reverse)
    save_pdb_cif_files(new_ent, pocket_name, output_pdb_file, output_cif_file)

def save_pocket_star(input_args):
    return save_pocket(*input_args)

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--pocket_names_file", type=str, required=True)
    parser.add_argument("--output_folder", type=str, required=True)
    parser.add_argument("--compound_lib", type=str, required=True)
    parser.add_argument("--num_threads", type=int, default=1)
    parser.add_argument("--max_protein_chains", type=int, default=10)
    parser.add_argument("--max_ligand_chains", type=int, default=10)
    parser.add_argument(
        "--cif_dir", type=str, default="/scicore/data/managed/PDB/latest/data/structures/divided/mmCIF/",
        help="Path to mmcif directory."
    )
    args = parser.parse_args()
    pocket_names = []
    with open(args.pocket_names_file) as f:
        for line in f:
            pocket_name = line.strip()
            pdb_id, biounit_id, protein_chains, ligand_chains, _ = pocket_name.split("__")
            protein_chains, ligand_chains = protein_chains.split("_"), ligand_chains.split("_")
            if len(protein_chains) > args.max_protein_chains or len(ligand_chains) > args.max_ligand_chains:
                logging.info(f"Too many chains for {pocket_name}")
                continue
            pocket_names.append(pocket_name)
    output_folder, cif_folder, compound_lib = args.output_folder, args.cif_dir, args.compound_lib
    Path(output_folder).mkdir(exist_ok=True)
    (Path(output_folder) / "fasta_files").mkdir(exist_ok=True)
    (Path(output_folder) / "pdb_files").mkdir(exist_ok=True)
    (Path(output_folder) / "cif_files").mkdir(exist_ok=True)
    (Path(output_folder) / "ligand_files").mkdir(exist_ok=True)
    input_args = [(pocket_name, cif_folder, output_folder, compound_lib) for pocket_name in tqdm(pocket_names)]
    num_threads = int(args.num_threads)
    
    with Pool(num_threads) as p:
        _ = list(tqdm(p.imap(save_pocket_star, input_args), total=len(input_args)))

if __name__ == "__main__":
    main()