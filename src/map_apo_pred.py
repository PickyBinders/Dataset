from pathlib import Path
from ost import io
from promod3.modelling import FSStructureServer, PentaMatch
import logging
from tqdm import tqdm
from collections import defaultdict
from create_dataset import load_cif_data, read_df

def get_afdb_structure(uniprot_id, fs_server):
    omf_obj = fs_server.GetOMF(uniprot_id)
    ent = omf_obj.GetAU()
    return ent    

def get_afdb_structure_from_idx(idx, fs_server):
    omf_obj = fs_server.GetOMFByIdx(idx)
    ent = omf_obj.GetAU()
    return ent    

def map_afdb(df_file, cif_data_dir, cif_dir, output_dir, fs_server, pentamatch_db_dir, use_pentamatch=False, overwrite=False, top_n=5):
    afdb_folder = Path(output_dir) / "afdb"
    afdb_folder.mkdir(exist_ok=True)
    pdb_chain_to_uniprot = load_cif_data(cif_data_dir, names_to_load=["uniprot_ids"])["uniprot_ids"]
    fs_server = FSStructureServer(fs_server)
    
    df = read_df(df_file)
    pdb_chains = defaultdict(set)
    uniprot_ids = set()
    for i, row in tqdm(df[["pocket_ID", "num_chains_in_pocket"]].drop_duplicates().iterrows()):
        if row["num_chains_in_pocket"] == 1:
            pdb_id, _, protein_chain, _, _ = row["pocket_ID"].split("__")
            uniprot_id = pdb_chain_to_uniprot.get(f"{pdb_id}_{protein_chain}", None)
            if uniprot_id is None:
                pdb_chains[pdb_id.lower()].add(protein_chain)
            else:
                uniprot_ids.add(uniprot_id.upper())

    for uniprot_id in tqdm(uniprot_ids):
        try:
            output_file = afdb_folder / f"{uniprot_id}.pdb"
            if output_file.exists() and not overwrite:
                continue
            afdb_entry = get_afdb_structure(uniprot_id, fs_server)
            io.SavePDB(afdb_entry, str(output_file))
        except Exception as e:
            logging.info(f"Couldn't get AFDB structure for {uniprot_id}: {e}")
    if not use_pentamatch:
        return
    
    pentamatch = PentaMatch(pentamatch_db_dir)
    afdb_hits = set()
    for pdb_id in tqdm(pdb_chains):
        cif_file = cif_dir / pdb_id[1:3] / f"{pdb_id}.cif.gz"
        if not cif_file.exists():
            continue
        _, seqres = io.LoadMMCIF(str(cif_file), seqres=True, info=False, remote=False)
        for s in seqres:
            if s.name in pdb_chains[pdb_id]:
                try:
                    afdb_hits |= set(pentamatch.TopN(top_n, s.string))
                except Exception as e:
                    logging.info(f"Couldn't get AFDB hits for {pdb_id}: {e}")
    
    for afdb_index in tqdm(afdb_hits):
        try:
            afdb_entry = get_afdb_structure_from_idx(afdb_index, fs_server)
            uniprot_id = afdb_entry.GetName().split()[0]
            output_file = afdb_folder / f"{uniprot_id}.pdb"
            if output_file.exists() and not overwrite:
                continue
            io.SavePDB(afdb_entry, str(output_file))
        except Exception as e:
            logging.info(f"Couldn't get AFDB structure for {uniprot_id}: {e}")


def map_apo(df_file, fs_lookup_file, output_dir):
    output_dir.mkdir(exist_ok=True)
    
    df = read_df(df_file)
    holo_chains = set(chain_key for pocket_id in df['pocket_ID'] for chain_key in [f"{pocket_id.split('__')[0].lower()}.cif.gz" + '_' + chain for chain in pocket_id.split('__')[2].split('_')])
    ligand_chains = set(chain_key for pocket_id in df['pocket_ID'] for chain_key in [f"{pocket_id.split('__')[0].lower()}.cif.gz" + '_' + chain for chain in pocket_id.split('__')[3].split('_')])
    oligo_states = dict(zip(df['PDB_ID'], df['oligo_state']))

    pdb_chains_all = set()
    with open(fs_lookup_file) as f:
        for line in tqdm(f):
            entry = line.strip().split('\t')[1]
            if "MODEL" in entry:
                continue
            pdb_chains_all.add(entry)
    
    apo_chains = pdb_chains_all - holo_chains - ligand_chains
    homomeric_chains = set()
    for pdb_id_chain in tqdm(apo_chains):
        pdb_id = pdb_id_chain.split(".cif.gz")[0].upper()
        if oligo_states.get(pdb_id, '').startswith("homo"):
            homomeric_chains.add(pdb_id_chain)
    apo_chains -= homomeric_chains
    
    with open(output_dir / "apo_chains.txt", "w") as f:
        for apo_chain in apo_chains:
            f.write(f"{apo_chain}\n")


    
def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--df_file", type=Path, required=True,
    )
    parser.add_argument(
        "--fs_lookup_file", type=Path, required=True,
    )
    parser.add_argument(
        "--cif_dir", type=Path, required=True,
    )
    parser.add_argument(
        "--cif_data_dir", type=Path, required=True,
    )
    parser.add_argument(
        "--output_dir", type=Path, required=True,
    )
    parser.add_argument(
        "--fs_server", type=str,
        help="""Path to the fs server containing structures from AFDB""",
    )
    parser.add_argument(
        "--pentamatch_db_dir", type=str,
        help="""Path to the pentamatch database containing sequences from AFDB""",
    )
    parser.add_argument(
        "--use_pentamatch", action="store_true",
    )
    parser.add_argument(
        "--overwrite", action="store_true",
    )
    args = parser.parse_args()

    map_afdb(args.df_file, args.cif_data_dir, args.cif_dir, args.output_dir, args.fs_server, args.pentamatch_db_dir, overwrite=args.overwrite, use_pentamatch=args.use_pentamatch)
    map_apo(args.df_file, args.fs_lookup_file, args.output_dir)

if __name__ == "__main__":
    main()