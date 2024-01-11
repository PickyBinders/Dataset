from pathlib import Path
from ost import io, mol
from promod3.modelling import FSStructureServer
import logging
import pandas as pnd
from tqdm import tqdm
import json
from collections import defaultdict
from create_dataset import read_df, load_cif_data, get_uniprot_to_chains

def get_apo_chain(pdb_id, chain_id, cif_dir):
    cif_file = Path(cif_dir) / pdb_id[1:3] / f"{pdb_id}.cif.gz"
    ent_full = io.LoadMMCIF(str(cif_file))
    ent =  mol.CreateEntityFromView(ent_full.Select(f"cname='{chain_id}'"), True)
    edi = ent.EditXCS(mol.BUFFERED_EDIT)
    for i, chain in enumerate(ent.GetChainList()):
        edi.RenameChain(chain, "A")
    edi.UpdateICS()
    return ent

def get_afdb_structure(uniprot_id, fs_server):
    omf_obj = fs_server.GetOMF(uniprot_id)
    ent = omf_obj.GetAU()
    return ent    

def map_apo_and_afdb(df_file, cif_data_dir, output_dir, cif_dir, fs_server, overwrite=False):
    pdb_chain_to_uniprot = load_cif_data(cif_data_dir, names_to_load=["uniprot_ids"])["uniprot_ids"]
    uniprot_to_pdb_chains = get_uniprot_to_chains(pdb_chain_to_uniprot)
    fs_server = FSStructureServer(fs_server)
    output_dir.mkdir(exist_ok=True)
    
    df = read_df(df_file)
    holo_chains = set(chain_key for pocket_id in df['pocket_ID'] for chain_key in [pocket_id.split('__')[0] + '_' + part for part in pocket_id.split('__')[2].split('_')])
    oligo_states = dict(zip(df['PDB_ID'], df['oligo_state']))

    uniprot_ids = set()
    for i, row in tqdm(df[["pocket_ID", "num_chains_in_pocket"]].drop_duplicates().iterrows()):
        if row["num_chains_in_pocket"] == 1:
            pdb_id, _, protein_chain, _, _ = row["pocket_ID"].split("__")
            uniprot_id = pdb_chain_to_uniprot.get(f"{pdb_id}_{protein_chain}", None)
            if uniprot_id is None:
                continue
            uniprot_ids.add(uniprot_id.upper())

    for uniprot_id in tqdm(uniprot_ids):
        uniprot_folder = output_dir / uniprot_id
        if uniprot_folder.exists() and not overwrite:
            continue
        uniprot_folder.mkdir(exist_ok=True)
        try:
            af2_ent = get_afdb_structure(uniprot_id, fs_server)
            io.SavePDB(af2_ent, str(uniprot_folder / f"{uniprot_id}_afdb.pdb"))
        except Exception as e:
            logging.info(f"Couldn't save AFDB structure for {uniprot_id}: {e}")
            continue
        for chain_key in uniprot_to_pdb_chains[uniprot_id]:
            if chain_key in holo_chains:
                continue
            pdb, chain = chain_key.split('_')
            if not oligo_states.get(pdb, '').startswith('homo'):
                try:
                    apo_ent = get_apo_chain(pdb.lower(), chain, cif_dir)
                    io.SavePDB(apo_ent, str(uniprot_folder / f"{uniprot_id}_{chain_key}.pdb"))
                except Exception as e:
                    logging.info(f"Couldn't save apo structure for {uniprot_id}_{chain_key}: {e}")
                    continue
    
def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--df_file", type=Path, required=True,
    )
    parser.add_argument(
        "--cif_data_dir", type=Path, required=True,
    )
    parser.add_argument(
        "--output_dir", type=Path, required=True,
    )
    parser.add_argument(
        "--fs_server", type=str, required=True,
        help="""Path to the fs server cointaining structures from AFDB""",
    )
    parser.add_argument(
        "--cif_dir", type=str, default="/scicore/data/managed/PDB/latest/data/structures/divided/mmCIF/",
        help="Path to mmcif directory."
    )
    parser.add_argument(
        "--overwrite", action="store_true",
    )
    args = parser.parse_args()

    map_apo_and_afdb(args.df_file, args.cif_data_dir, args.output_dir, args.cif_dir, args.fs_server, overwrite=args.overwrite)

if __name__ == "__main__":
    main()