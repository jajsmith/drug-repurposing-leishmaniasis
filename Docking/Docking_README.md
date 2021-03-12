# Dependencies

### Software (Linux)
* `Open Babel`
* `AutoDock Vina` 

### Python packages
* `OpenBabel + Pybel`
* `BioPython`
* `Numpy`
* `Pandas`

# Instructions

## Step 1 - Generate protein files

```bash
time python google_cloud_proteins/extract_google_pdb_files.py

time python rcsb_proteins/extract_rcsb_pdb_files.py
```

After running the two commands above, you will see the following folders for each type of protein file you will later use for docking with AutoDock Vina:

* `PDBs` will contain the raw (not docking ready) PDB files for each Leishmania protein found in either INDABA's Google Cloud folder or in the RCSB database (as of February 1 2021) 
  
* `single-chain-PDBs` will contain docking-ready PDB files for each Leishmania protein from the aforementioned sources. If a raw PDB file in `PDBs` contains multiple units of a protein plus a ligand and water atoms, then the corresponding single chain PDB file will contain only one chain (usually chain A) of the associated raw PDB file

* `Sequences` will contain a single FASTA file with the amino acid sequences of the chains inside `single-chain-PDBs`

* `PDBQTs` will contain the PDBQT files of each of the docking-ready chains in `single-chain-PDBs`. These .pdbqt files are the actual input to AutoDock vina when docking is run

* `CONFIGs` will contain AutoDock configuration files for each of the proteins in `single-chain-PDBs`. By default, the coordinates of the search box are centered at the center of the protein atoms, and the size is calculated such that the box encloses the entire protein. In other words, these config files work for "unconstrained search" when running docking with AutoDock vina

## Step 2 - Generate ligand files

```bash
time python indaba_drugs/extract_indaba_drug_files.py
```

After running the above command, you will see the following folder

* `Ligands` will contain .pdbqt files of all drugs included in the Indaba challenge (in-trials, Drug Central, and DrugBank). These will be used as inputs on the docking runs with AutoDock vina

## Step 3 - Run docking screening
TODO