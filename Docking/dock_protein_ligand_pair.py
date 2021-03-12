import openbabel as ob
from openbabel import pybel
import pandas as pd
import argparse
import pyrosetta
import os
import sys
import subprocess

"""
@Description
This script docks a protein-ligand pair given a Protein ID, a Drug ID, 
and the location of the AutoDock vina executable file. Results will be saved
inside the Docking-Results folder with the name <ProteinID>_<DrugID>.pdbqt
and the log of the docking run as <ProteinID>_<DrugID>.log

@Example
python dock_protein_ligand_pair.py -p 6RXC -l ZINC000029747110 -v ~/autodock_vina/bin/vina

@Note
Please ensure that you pass the location of your vina executable file to the "v" argument

@Instructions
You should first run the Python scripts to generate docking-ready PDBQT files for both proteins and drugs (see Docking_README.md for instructions)

csv_dir = args.i
out_filename = args.o
"""

# Configure argparser and load arguments
parser = argparse.ArgumentParser(description="Dock a protein-ligand pair with AutoDock vina")
parser.add_argument("-p", help="Protein identifier (e.g. 6AZ1)", default="6RXC")
parser.add_argument("-l", help="Ligand or drug identifier (e.g. DB11989)", default="DB00336")
parser.add_argument("-v", help="Path to AutoDock Vina executable", default="~/Downloads/autodock_vina_1_1_2_linux_x86/bin/vina")
args = parser.parse_args()

""" TODO Finish this!
VINA_PATH = "~/autodock_vina_1_1_2_linux_x86/bin/vina" # NOTE change this to the Vina path in your system
PDB_DIR = "../Proteins/PDB"
PDBQT_DIR = "./Proteins/PDBQT"
CONFIG_DIR = "./Configs"
LIGANDS_DIR = "./Ligands"
PROTEINS_DIR = "./Proteins"
OUTPUT_DIR = "vina-results"

if not os.path.exists("./Docking-Results"):
    os.mkdir("./Docking-Results")

# Get ligands and proteins
f = open(f"{LIGANDS_DIR}/ligands.txt", "r")
ligands = f.read().splitlines()
f.close()

f = open(f"{PROTEINS_DIR}/proteins.txt", "r")
proteins = f.read().splitlines()
f.close()

results_ligands = []
results_proteins = []
results_vinascores = []

for ligand in ligands:
    for protein in proteins:
        results_ligands.append(ligand)
        results_proteins.append(protein)

        # Get receptor PDBQT file
        receptor_file = f"{PDBQT_DIR}/{protein}sc.clean.pdbqt"
        
        # Get ligand PDBQT file

        ligand_file = f"{LIGANDS_DIR}/{ligand}.pdbqt"

        # Get Vina configuration file
        vina_config_file = f"{CONFIG_DIR}/{protein}sc_vina.config"

        # Get Vina output files
        log_file = f"{OUTPUT_DIR}/{protein}-{ligand}_vina_log.txt"
        out_file = f"{OUTPUT_DIR}/{protein}-{ligand}_vina_results.pdbqt"
        command = f"{VINA_PATH} --receptor {receptor_file} --ligand {ligand_file} --config {vina_config_file} --out {out_file} --log {log_file}"
        os.system(command)

        # Add top score to file report
        try:
            f = open(log_file, 'r')
            lines = f.read().splitlines()
            f.close()
            top_result_index = 0
            for n, line in enumerate(lines):
                if "-----+" in line:
                    top_result_index = n + 1
                    break
            line = lines[top_result_index]
            start = line.index('-')
            end = line.index(' ', start)
            score = line[start:end]
            results_vinascores.append(score)
            
        except:
            results_vinascores.append("0")
            continue

out_df = pd.DataFrame()
out_df['protein'] = results_proteins
out_df['ligand'] = results_ligands
out_df['vina-score'] = results_vinascores
out_df.to_csv("screening-results.csv", index=False)

print("FINISHED WITH SUCCESS!")"""