import json
from ost import conop, mol, io
from pathlib import Path
from tqdm import tqdm
import logging
from create_dataset import read_df

# Define available names for protein and ligand chains
PROTEIN_CHAINS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
LIGAND_CHAINS = PROTEIN_CHAINS.lower() + '0123456789'

def rename_chains(ent, ent_selected, protein_chains, ligand_chains):
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
        else:
            logging.info(f"Chain {original_name} not in protein or ligand chains")
            return None, None
        edi.RenameChain(chain, final_name)
        edi.SetChainDescription(chain, original_chain.description)
        edi.SetChainType(chain, original_chain.type)
        name_mapping[original_name] = final_name
    edi.UpdateICS()
    return name_mapping, ent_selected

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
    return biounit, seqres, info

def save_ligands(ent, ligand_chains, output_prefix):
    for ligand_chain in ligand_chains:
        ligand = ent.Select(f"chain='{ligand_chain}'")
        try:
            io.SaveEntity(ligand, str(output_prefix / f"{ligand_chain}.sdf"))
        except Exception as e:
            logging.error(f"Could not save ligand {ligand_chain} in {output_prefix}: {e}")
            return

def save_sequences(ent, output_fasta_file, chain_to_sequence):
     with open(output_fasta_file, "w") as f:
        for x in ent.chains:
            if x.name in chain_to_sequence:
                f.write(f">{x.name}\n")
                f.write(chain_to_sequence[x.name] + "\n")

def save_pdb_file(ent_renamed, output_pdb_file, name_mapping, output_mapping_file):
    try:
        io.SavePDB(ent_renamed, str(output_pdb_file))
        with open(output_mapping_file, "w") as f:
            json.dump(name_mapping, f)
    except Exception as e:
        logging.error(f"Could not save {output_pdb_file}: {e}")

def save_cif_file(ent, output_cif_file, info, pocket_name):
    try:
        lib = conop.GetDefaultLib()
        entity_info = io.MMCifWriterEntityList()
        entity_ids = set(info.GetMMCifEntityIdTr(ch.name.split(".")[-1]) for ch in ent.chains)
        for entity_id in info.GetEntityIdsOfType("polymer"):
            if entity_id not in entity_ids:
                continue
            # Get entity description from info object
            entity_desc = info.GetEntityDesc(entity_id)
            # interface of entity_desc is similar to MMCifWriterEntity
            entity_poly_type = entity_desc.entity_poly_type
            mon_ids = entity_desc.mon_ids
            e = io.MMCifWriterEntity.FromPolymer(entity_poly_type, mon_ids, lib)
            entity_info.append(e)
            # search all chains assigned to the entity we just added
            for ch in ent.chains:
                if info.GetMMCifEntityIdTr(ch.name.split(".")[-1]) == entity_id:
                    entity_info[-1].asym_ids.append(ch.name)
            # deal with heterogeneities
            for a,b in zip(entity_desc.hetero_num, entity_desc.hetero_ids):
                entity_info[-1].AddHet(a,b)
        writer = io.MMCifWriter()
        writer.SetStructure(ent, lib, entity_info=entity_info)
        writer.Write(pocket_name, str(output_cif_file))    
    except Exception as e:
        logging.error(f"Could not save {output_cif_file}: {e}")

def save_pocket(pocket_name, cif_folder, output_folder, overwrite=False):
    if len(pocket_name) > 200:
        logging.error(f"{pocket_name} too long for file creation. Skipping.")
        return
    pdb_id, biounit_id, protein_chains, ligand_chains, _ = pocket_name.split("__")
    protein_chains, ligand_chains = protein_chains.split("_"), ligand_chains.split("_")
    cif_file = Path(cif_folder) / f"{pdb_id.lower()[1:3]}" / f"{pdb_id.lower()}.cif.gz"
    if not cif_file.exists():
        logging.error(f"{cif_file} does not exist")
        return
    output_folder_pocket = Path(output_folder) / pocket_name
    output_folder_pocket.mkdir(exist_ok=True)
    output_pdb_file = output_folder_pocket / f"system.pdb"
    output_name_mapping_file = output_folder_pocket / f"chain_mapping.json"
    output_cif_file = output_folder_pocket / f"system.cif"
    output_fasta_file = output_folder_pocket / f"sequences.fasta"
    output_ligand_prefix = output_folder_pocket / "ligand_files"
    output_ligand_prefix.mkdir(exist_ok=True)
    if output_cif_file.exists() and not overwrite:
        return
    biounit, seqres, info = load_biounit_seqres(cif_file, biounit_id)
    if biounit is None or seqres is None:
        return
    chain_to_sequence = {s.name: s.string for s in seqres}
    chain_to_sequence = {x: chain_to_sequence[x.split(".")[-1]] for x in protein_chains if x.split(".")[-1] in chain_to_sequence}
    ent_selected = mol.CreateEntityFromView(biounit.Select(" or ".join(f"chain='{c}'" for c in protein_chains + ligand_chains)), True)
    save_ligands(ent_selected, ligand_chains, output_ligand_prefix)
    save_sequences(ent_selected, output_fasta_file, chain_to_sequence)
    save_cif_file(ent_selected, output_cif_file, info, pocket_name)
    if len(protein_chains) > len(PROTEIN_CHAINS) or len(ligand_chains) > len(LIGAND_CHAINS):
        logging.error(f"Too many chains to convert to PDB for {pocket_name}")
        return
    name_mapping, ent_renamed = rename_chains(biounit, ent_selected.Copy(), protein_chains, ligand_chains)
    save_pdb_file(ent_renamed, output_pdb_file, name_mapping, output_name_mapping_file)

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--df_file", type=str, required=True)
    parser.add_argument("--output_folder", type=str, required=True)
    parser.add_argument("--midfix", type=str, default=None)
    parser.add_argument(
        "--cif_dir", type=str, default="/scicore/data/managed/PDB/latest/data/structures/divided/mmCIF/",
        help="Path to mmcif directory."
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    pocket_names = set(read_df(args.df_file)["system_ID"])
    if args.midfix is not None:
        pocket_names = set(x for x in pocket_names if x[1:3].lower() == args.midfix)
    output_folder, cif_folder = args.output_folder, args.cif_dir
    Path(output_folder).mkdir(exist_ok=True)
    if len(pocket_names) == 0:
        return
    for pocket_name in tqdm(pocket_names):
        save_pocket(pocket_name, cif_folder, output_folder, args.overwrite)

if __name__ == "__main__":
    main()