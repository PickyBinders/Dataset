from ost import io
from ost import mol

from pathlib import Path
import json
import logging
from tqdm import tqdm
    

def extract_pockets(cif_dir, out_file, threshold=8):
    """
    Extracts binding site residues and center of mass for each ligand in each structure in the given directory.
    The binding site residues are defined as all residues within a given distance threshold from the ligand.
    Key for the returned dictionary is the PDB ID, ligand name and ligand MMCIF chain.
    """
    cifs_dir = Path(cif_dir)
    data = {}

    for cif_file in tqdm(cifs_dir.iterdir()):
        plc, info = io.LoadMMCIF(str(cif_file), info=True)
        entry_id = cif_file.stem.split(".")[0].upper()
        resnum_mapping = dict()
        protein_chains = [chain for chain in plc.chains if chain.type == mol.CHAINTYPE_POLY_PEPTIDE_L]
        for protein_chain in protein_chains:
            resnum_mapping[protein_chain.name] = dict([(k, v) for v, k in enumerate([r.number for r in protein_chain.residues])])
        num = 0
        ligand_chains = [chain for chain in plc.chains if chain.type == mol.CHAINTYPE_NON_POLY]
        for ligand_chain in ligand_chains:
            for ligand in ligand_chain.residues:
                binding_site = plc.Select(f"{threshold} <> [cname={ligand_chain.name} and rnum={ligand.number.num}] and protein=True")
                try:
                    binding_site = [{"chain": r.chain.name, 
                           "one_letter_code": r.one_letter_code, 
                           "residue_index": resnum_mapping[r.chain.name][r.number.num], 
                           "residue_number": r.number.num} for r in binding_site.residues]
                    num += 1
                except Exception as e:
                    binding_site = None
                    logging.exception(f'Extraction of binding site residues failed for {cif_file}')
                if binding_site is not None and len(binding_site) > 0:
                    data[f'{entry_id}_{ligand.name}_{ligand_chain.name}'] = {"binding_site": binding_site, "center_of_mass": list(ligand.GetCenterOfMass())}
        print(f"Extracted {num} binding sites for {entry_id}")
    with open(out_file, "w") as f:
        json.dump(data, f)


def main():
    from sys import argv
    cif_dir, out_file = argv[1:]
    extract_pockets(cif_dir, out_file)

if __name__ == "__main__":
    main()