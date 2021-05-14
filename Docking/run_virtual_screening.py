import argparse
import os
import sys
import subprocess

"""
@Description
This script screens all protein-molecule pairs in a given comma-separated file.
The input .csv file should not have a header and each row should only have
the <protein id> followed by the <ligand/drug id> separated by a comma. See template below:

6RXC,ZINC000029747110
6AZ1,DB06243
2P1E,DC4
A0A088RHB8,DB11820

@Example
python run_virtual_screening -i input_pairs.csv -v ~/Downloads/autodock_vina_1_1_2_linux_x86/bin/vina
"""

# Configure argparser and load arguments
parser = argparse.ArgumentParser(description="Run a virtual screening given a list of protein-ligand pairs")
parser.add_argument("-i", help="Input .csv file with protein-ligand pairs")
parser.add_argument("-v", help="Path to AutoDock Vina executable", default="~/Downloads/autodock_vina_1_1_2_linux_x86/bin/vina")
args = parser.parse_args()

INPUT_FILE = args.i
VINA_PATH = args.v

# Get protein-ligand pairs
f = open(INPUT_FILE, "r")
pairs = f.read().splitlines()
f.close()

for pair in pairs:
    protein, ligand = pair.split(',')
    command = f"python dock_protein_ligand_pair.py -p {protein} -l {ligand} -v {VINA_PATH}"
    os.system(command)

print("FINISHED WITH SUCCESS!")