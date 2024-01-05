from ost import io
from ost import mol

from pathlib import Path
import json
import logging
from tqdm import tqdm

def extract_cif_data(cif_dir, out_file, threshold=8):
    """
    Extracts binding site residues and center of mass for each ligand in each structure in the given directory.
    The binding site residues are defined as all residues within a given distance threshold from the ligand.
    Key for the returned dictionary is the PDB ID, ligand name and ligand MMCIF chain.
    """
    cifs_dir = Path(cif_dir)
    data = dict(dates=dict(), pockets=dict())
    for cif_file in tqdm(cifs_dir.iterdir()):
        try:
            plc, info = io.LoadMMCIF(str(cif_file), info=True)
        except Exception as e:
            logging.error(f"Could not load {cif_file}: {e}")
            continue
        entry_id = info.struct_details.entry_id.upper()
        data["dates"][entry_id] = info.revisions.GetDateOriginal()
        protein_chains = [chain for chain in plc.chains if chain.type == mol.CHAINTYPE_POLY_PEPTIDE_L]
        if len(protein_chains) == 0:
            continue
        resnum_mapping = dict()
        for protein_chain in protein_chains:
            resnum_mapping[protein_chain.name] = dict([(k, v) for v, k in enumerate([r.number for r in protein_chain.residues])])
        ligand_chains = [chain for chain in plc.chains if chain.type == mol.CHAINTYPE_NON_POLY]
        if len(ligand_chains) == 0:
            continue
        num = 0
        ligand_selections = {}
        for ligand_chain in ligand_chains:
            for ligand in ligand_chain.residues:
                selection = plc.Select(f"{threshold} <> [cname={ligand_chain.name} and rnum={ligand.number.num}] and protein=True")
                try:
                    pocket = [{"chain": r.chain.name, 
                            "one_letter_code": r.one_letter_code, 
                            "residue_index": resnum_mapping[r.chain.name][r.number.num], 
                            "residue_number": r.number.num} for r in selection.residues]
                except Exception as e:
                    pocket = None
                if pocket is not None and len(pocket) > 0:
                    ligand_selections[(ligand_chain.name, ligand.name)] = pocket
        for biounit in info.biounits:
            biounit_chains = set(str(c) for c in biounit.GetChainList())
            for ligand_chain in ligand_chains:
                if ligand_chain.name not in biounit_chains:
                    continue
                for ligand in ligand_chain.residues:
                    if (ligand_chain.name, ligand.name) not in ligand_selections:
                        continue
                    binding_site = [r for r in ligand_selections[(ligand_chain.name, ligand.name)] if r["chain"] in biounit_chains]
                    if len(binding_site) > 0:
                        num += 1
                        data["pockets"][f'{entry_id}__{biounit.id}__{ligand_chain.name}'] = {"pocket_residues": binding_site, "center_of_mass": list(ligand.GetCenterOfMass())}
        if num > 0:
            print(f"Extracted {num} binding sites for {entry_id}")
    with open(out_file, "w") as f:
        json.dump(data, f)


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("cif_dir", help="Directory containing the CIF files")
    parser.add_argument("out_file", help="Output file")
    parser.add_argument("--threshold", help="Distance threshold for binding site residues", default=8)
    args = parser.parse_args()
    extract_cif_data(args.cif_dir, args.out_file, args.threshold)

if __name__ == "__main__":
    main()