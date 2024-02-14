# Protein-Ligand Complex Prediction Dataset

## Available Data
- `data_new.tsv`: The main dataset, containing all ligand systems and their annotations (previous version in `data.tsv` and `data_clean.tsv`)
- `column_info.csv`: A description of each column in `data_new.tsv`
- `component_mapping_new.tsv`: The connected component IDs for each system
- `system_files_new`: The folder containing all files related to each system in corresponding subfolders. Each system subfolder has, `system.cif`, `system.fasta`, a folder called `ligand_files` with SDF files of each ligand, `system.pdb`, and `chain_mapping.json` describing the mapping between the original chain and renamed chain in the .pdb file.
- `scores`, `components`, `component_info.csv`: Results of similarity searches and graph clustering
- `validation`: Results of extracting validation report information
- `cif_data`: Results of extracting data from CIF files
- `annotations`: Examples of SMTL annotation files with PLIP data
- `cross_docking`: Scores linking system_IDs to each other for cross-docking
- `apo_scores`: Scores linking apo chains to system_IDs
- `afdb_scores`: Scores linking AFDB structures to system_IDs

The Python scripts referenced below are in https://github.com/PickyBinders/Dataset

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
- The final dataset consists of all systems with at least one non-ion and non-artifact ligand having <10 protein chains and <10 ligand chains
See `column_info.csv` for a description of each resulting column in the dataframe

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

This resource is developed by the Computational Structural Biology Group at the SIB Swiss Institute of Bioinformatics and the Biozentrum of the University of Basel. All copyrightable parts of this resource are available under a Creative Commons Attribution-ShareAlike 4.0 International License (see `LICENSE`)