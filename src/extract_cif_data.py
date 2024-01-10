from ost import io
from ost import mol

from pathlib import Path
import json
import logging
from tqdm import tqdm

def extract_cif_data(cif_dir, out_file, ignore, threshold=8):
    """
    Extracts dates for each PDB ID and binding site residues and center of mass for each ligand in each structure in the given directory.
    The binding site residues are defined as all residues within a given distance threshold from the ligand.
    Binding sites are defined for each ligand chain in each biounit.
    Key for the returned dictionary is the PDB ID, biounit and ligand MMCIF chain.
    """
    cifs_dir = Path(cif_dir)
    data = dict(dates=dict(), pockets=dict())
    for cif_file in tqdm(cifs_dir.iterdir()):
        try:
            plc, info = io.LoadMMCIF(str(cif_file), info=True)
        except Exception as e:
            logging.error(f"PDBError: Could not load {cif_file}: {e}")
            continue
        entry_id = info.struct_details.entry_id.upper()
        if entry_id in ignore:
            continue
        data["dates"][entry_id] = info.revisions.GetDateOriginal()
        num = 0
        for b in info.biounits:
            try:
                biounit = mol.alg.CreateBU(plc, b)
            except Exception as e:
                logging.error(f"BioUnitError: Could not create biounit for {cif_file}: {e}")
                continue
            protein_chains = [chain for chain in biounit.chains if chain.type == mol.CHAINTYPE_POLY_PEPTIDE_L]
            if len(protein_chains) == 0:
                continue
            ligand_chains = [chain for chain in biounit.chains if chain.type == mol.CHAINTYPE_NON_POLY]
            if len(ligand_chains) == 0:
                continue
            resnum_mapping = dict()
            for protein_chain in protein_chains:
                resnum_mapping[protein_chain.name] = dict([(k, v) for v, k in enumerate([r.number for r in protein_chain.residues])])
            for ligand_chain in ligand_chains:
                for ligand in ligand_chain.residues:
                    selection = biounit.Select(f"{threshold} <> [cname='{ligand.chain.name}' and rnum={ligand.number.num}] and protein=True")
                    try:
                        pocket = [{"chain": r.chain.name.split(".")[-1], 
                                "one_letter_code": r.one_letter_code, 
                                "residue_index": resnum_mapping[r.chain.name][r.number.num], 
                                "residue_number": r.number.num} for r in selection.residues]
                    except Exception as e:
                        pocket = None
                    if pocket is not None and len(pocket) > 0:
                        num += 1
                        ligand_chain_name = ligand.chain.name.split(".")[-1]
                        data["pockets"][f'{entry_id}__{b.id}__{ligand_chain_name}'] = {"pocket_residues": pocket, "center_of_mass": list(ligand.GetCenterOfMass())}
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
    parser.add_argument("--ignore", help="PDB IDs to ignore", nargs="+", default=[]) # to ignore some pathological cases
    args = parser.parse_args()
    extract_cif_data(args.cif_dir, args.out_file, args.ignore, args.threshold)

if __name__ == "__main__":
    main()