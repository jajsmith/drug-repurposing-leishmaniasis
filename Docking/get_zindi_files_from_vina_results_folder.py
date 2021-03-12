"""
Arguments

1 - <vina_results_dir> Directory with .pdbqt files resulting from Vina
2 - <protein_pdb_dir> Input folder with PDBs of proteins
6 - <output_dir> Output Directory to save the ZINDI-ready .pdb files
"""

import os
import sys
import pandas as pd

vina_dir = sys.argv[1]
protein_pdb_dir = sys.argv[2]
output_dir = sys.argv[3]

# Get ligand names
vina_files = os.listdir(vina_dir)
vina_files = [file for file in vina_files if '.pdbqt' in file]
proteins = [file.split("-")[0] for file in vina_files]

PDB_DIR = protein_pdb_dir

for i in range(len(vina_files)):
    pdbqt_file = vina_files[i]
    pdb_file = f"{PDB_DIR}/{proteins[i]}sc.clean.pdb"
    print(f"Converting file {pdbqt_file}")
    command = f"python prepare_zindi_file.py {vina_dir} {pdbqt_file} {output_dir} {pdb_file}"
    os.system(command)

print("Finished with success!")
