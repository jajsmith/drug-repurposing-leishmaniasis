import openbabel as ob
from openbabel import pybel
import pubchempy as pcp
import pandas as pd
import argparse
import pyrosetta
import os
import sys
import subprocess

"""
@Description
This script screens all protein-molecule pairs given

# Configure argparser and load arguments
parser = argparse.ArgumentParser(description="Extract data from MarketScan tables and save it as a Python dictionary")
parser.add_argument("-i", help="Directory with all input .csv or .tsv MarketScan tables", default="/mnt/rstor/SOM_CCCC_XXB16/pranteshDATA/Jahir/CSV_data_ICI_patients/")
parser.add_argument("-o", help="Pickle file to save results in", default="/mnt/rstor/SOM_CCCC_XXB16/pranteshDATA/Jahir/ici_patients_manuscript.pkl")
args = parser.parse_args()

csv_dir = args.i
out_filename = args.o
"""

VINA_PATH = "/home/jahir/Desktop/RESEARCH_PROJECTS/LEISHMANIASIS/autodock_vina_1_1_2_linux_x86/bin/vina"
PDB_DIR = "../Proteins/PDB"
PDBQT_DIR = "./Proteins/PDBQT"
CONFIG_DIR = "./Configs"
LIGANDS_DIR = "./Ligands"
PROTEINS_DIR = "./Proteins"
OUTPUT_DIR = "vina-results"

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

print("FINISHED WITH SUCCESS!")