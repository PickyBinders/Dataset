# Protein-Ligand Complex Prediction Dataset
This is a reference implementation of protein-ligand dataset curation, annotation and similarity calculation. 
This effort is now continued in PLINDER (https://github.com/aivant/plinder/, to be open-sourced shortly).

## Data Sources
1. Protein Data Bank (PDB) MMCIF files (https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF)
See `extract_cif_data.py`. 
2. PDB XML validation reports (https://files.wwpdb.org/pub/pdb/validation_reports/)
See `extract_validation_data.py`.
3. Swiss Model Template Library (e.g https://swissmodel.expasy.org/templates/5KKH)
See `extract_plip_data.py`.
4. PDB Chemical Components Dictionary (https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz)
5. PDB seqres (https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt)
6. AlphaFold Database (https://alphafold.ebi.ac.uk/)
7. UniRef50 (https://www.ebi.ac.uk/uniprot/download-center)

## Merging and annotating systems
- All pockets and associated info from the 3 sources are merged into a single dataframe, using the entry_pdb_id, entry_biounit and ligand_chain as the key and annotated (see `create_dataset.py`). Each row of the resulting dataframe is referred to as a ligand system - identified by its (entry_pdb_id, entry_biounit, ligand_chain, ligand_interacting_protein_chains)
- ligand systems are merged into systems if they share the same entry_pdb_id and entry_biounit, and the ligands are within 4A of each other.
- Subsets are created by
    - filtering out systems with only ions
    - filtering out systems with only ions and/or artifacts
- The final dataset consists of all systems with at least one non-ion and non-artifact ligand having <10 protein chains and <10 ligand chains.

## System clustering

A Foldseek search is run on all PDB chains participating in a system (parameters in `run_foldseek.sh`), and the alignments and per-residue lDDT scores are saved. Protein, pocket, pli and ligand-level scores between systems are extracted as in `extract_scores.py`. For each individual score, a graph is created with nodes as systems and edges between systems with a score above a given threshold (see `graph_clustering.py`). Connected components are extracted from these graphs for each system.

## Pairing structures
Systems are paired to cross-docking cases, apo chains and predicted chains as described below. For all three similarity scores between the system and the linked structures are available at the protein and pocket level.

### Cross-docking
For each system, all systems from a different PDB ID which have (`protein_fident_qcov_foldseek_weighted_sum` or `protein_fident_qcov_mmseqs_weighted_sum` > 98%) and (`pocket_fident_foldseek_weighted_sum` or `pocket_fident_mmseqs_weighted_sum` > 98%) are collected as candidates for cross-docking.

### Apo structures
All PDB protein chains not participating in a pocket containing a ligand with >=5 rotatable bonds are considered as apo chains. Each system consisting of a single interacting protien chain is linked where possible to apo chains which have `protein_fident_qcov_mmseqs_weighted_sum` >98%. 

### Predicted structures
An Mmseqs search is run for all system chains against all UniRef50 members with an AFDB structure of the UniProt IDs associated with each holo chain with `--min-seq-id 0.98 -c 0.9`. Each system consisting of a single interacting protein chain is linked where possible to predicted structures which have `protein_fident_qcov_mmseqs_weighted_sum` >98%.

## Authors
Janani Durairaj
Xavier Robin

This resource is developed by the Computational Structural Biology Group at the SIB Swiss Institute of Bioinformatics and the Biozentrum of the University of Basel. 
This work was supported by the LIGATE project. This project has received funding from the European High-Performance Computing Joint Undertaking (JU) under grant agreement No 956137. The JU receives support from the European Union’s Horizon 2020 research and innovation programme and Italy, Sweden, Austria, the Czech Republic, Switzerland.

All copyrightable parts of this resource are available under an Apache 2.0 License (see `LICENSE`)
