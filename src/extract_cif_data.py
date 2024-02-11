from collections import defaultdict
from ost import io, mol

from pathlib import Path
import json
import logging
from tqdm import tqdm

def extract_cif_data(cif_dir, out_file, ignore, threshold=6):
    """
    Extracts dates for each PDB ID and binding site residues and center of mass for each ligand in each structure in the given directory.
    The binding site residues are defined as all residues within a given distance threshold from the ligand.
    Binding sites are defined for each ligand chain in each biounit.
    Key for the returned dictionary is the PDB ID, biounit and ligand MMCIF chain.
    """
    cifs_dir = Path(cif_dir)
    data = dict(dates=dict(), pockets=dict(), uniprot_ids=dict(), 
                chain_mapping=dict(), lengths=dict(), 
                entity_mapping=dict(), instance_starts=dict())
    for cif_file in tqdm(cifs_dir.iterdir()):
        try:
            plc, seqres, info = io.LoadMMCIF(str(cif_file), info=True, seqres=True)
        except Exception as e:
            logging.error(f"PDBError: Could not load {cif_file}: {e}")
            continue
        entry_id = info.struct_details.entry_id.upper()
        if entry_id in ignore:
            continue
        data["dates"][entry_id] = info.revisions.GetDateOriginal()
        data["lengths"][entry_id] = {x.name: len(x.string) for x in seqres}
        author_chain_mapping = {chain.GetStringProp("pdb_auth_chain_name"): chain.name for chain in plc.chains if chain.type == mol.CHAINTYPE_POLY_PEPTIDE_L}
        if len(author_chain_mapping):
            data["chain_mapping"][entry_id] = author_chain_mapping
        data["entity_mapping"][entry_id] = {c.GetName(): info.GetMMCifEntityIdTr(c.GetName()) for c in plc.chains}
        for struct_ref in info.struct_refs:
            if struct_ref.db_name == "UNP":
                for aligned_seq in struct_ref.aligned_seqs:
                    chain_name = author_chain_mapping.get(aligned_seq.chain_name, None)
                    if chain_name is not None:
                        data["uniprot_ids"][f"{entry_id}_{chain_name}"] = struct_ref.db_access
        num = 0
        for b in info.biounits:
            try:
                biounit = mol.alg.CreateBU(plc, b)
            except Exception as e:
                logging.error(f"BioUnitError: Could not create biounit for {cif_file}: {e}")
                continue
            instance_starts_biounit = defaultdict(list)
            for chain in biounit.chains:
                instance_starts_biounit[chain.name.split(".")[-1]].append(int(chain.name.split(".")[0]))
            data["instance_starts"][f"{entry_id}__{b.id}"] = {k: min(v) for k, v in instance_starts_biounit.items() if min(v) > 1}
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
                        pocket = [{"chain": r.chain.name, 
                                "one_letter_code": r.one_letter_code, 
                                "residue_index": resnum_mapping[r.chain.name][r.number.num], 
                                "residue_number": r.number.num} for r in selection.residues]
                    except Exception as e:
                        pocket = None
                    if pocket is not None and len(pocket) > 0:
                        num += 1
                        data["pockets"][f'{entry_id}__{b.id}__{ligand.chain.name}'] = {"pocket_residues": pocket, "center_of_mass": list(ligand.GetCenterOfMass())}
        if num > 0:
            print(f"Extracted {num} binding sites for {entry_id}")
    with open(out_file, "w") as f:
        json.dump(data, f)


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("cif_dir", help="Directory containing the CIF files")
    parser.add_argument("out_file", help="Output file")
    parser.add_argument("--threshold", help="Distance threshold for binding site residues", default=6)
    parser.add_argument("--ignore", help="PDB IDs to ignore", nargs="+", default=[]) # to ignore some pathological cases
    args = parser.parse_args()
    extract_cif_data(args.cif_dir, args.out_file, args.ignore, args.threshold)

if __name__ == "__main__":
    main()