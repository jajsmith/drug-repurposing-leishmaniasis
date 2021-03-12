# Dependencies

### OS
This pipeline has been tested on Linux Ubuntu 18.04

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

* `LIGANDs` will contain .pdbqt files of all drugs included in the Indaba challenge (in-trials, Drug Central, and DrugBank). These will be used as inputs on the docking runs with AutoDock vina

## Step 3 - Download AutoDock Vina

```bash
wget http://vina.scripps.edu/download/autodock_vina_1_1_2_linux_x86.tgz -O ~/Downloads/autodock_vina_1_1_2_linux_x86.tgz

tar -xvzf ~/Downloads/autodock_vina_1_1_2_linux_x86.tgz
```

## Step 4a - Run docking for a protein-ligand pair

You should now be able to dock any protein-ligand pair from the molecule files generated above. For example, you could try to dock the structure 6RXC with the in-trial drug ZINC000029747110 as follows:

```bash
python dock_protein_ligand_pair.py -p 6RXC -l ZINC000029747110
```

The .pdbqt containing the ligand poses found by vina will be stored in a folder called `Docking-Results`
Note that you can use the argument `-v` to specify a different location of the Vina binary executable if you decided to install AutoDock vina in a location other than `~/Downloads`

The top score from all poses found by vina will be recorded in the file called `Docking-Score-Report.csv`. Note that the size of this file will grow as you continue to run docking pairs

## Step 4b - Run a virtual screening with several protein-ligand pairs

If you want to dock several protein-drug pairs, you can do so by creating a comma-separated file of protein-ligand pairs (with no header, see `example_pairs.csv` for an example) and pass it as an argument to the `run_virtual_screening.py` script.

```bash
python run_virtual_screening.py -i example_pairs.csv
```

The .pdbqt containing the ligand poses found by vina will be stored in a folder called `Docking-Results`

Also, the top score from all poses found by vina for each pair will be recorded in the file called `Docking-Score-Report.csv`. Note that the size of this file will grow as you continue to run docking pairs

## Step 5 - Convert docking results to a Zindi-ready submission PDB file

Once you have docked the pairs you selected, you can create Zindi-ready submission files with 

```bash
python get_zindi_files_from_vina_results_folder.py
```

The submission-ready files will be stored in `Zindi-Ready-Files` in .pdb format. These can be uploaded to Zindi's competition portal for scoring