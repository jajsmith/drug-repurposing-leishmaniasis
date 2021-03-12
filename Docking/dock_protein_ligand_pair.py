import argparse
import os

"""

@Description
This script docks a protein-ligand pair given a Protein ID, a Drug ID, 
and the location of the AutoDock vina executable file. Results will be saved
inside the Docking-Results folder with the name <ProteinID>_<DrugID>.pdbqt
and the log of the docking run as <ProteinID>_<DrugID>.log

@Example
python dock_protein_ligand_pair.py -p 6RXC -l ZINC000029747110 -v ~/Downloads/autodock_vina_1_1_2_linux_x86/bin/vina

@Note
Please ensure that you pass the location of your vina executable file to the "v" argument

@Instructions
You should first run the Python scripts to generate docking-ready PDBQT files for both proteins and drugs (see Docking_README.md for instructions)

"""

# Configure argparser and load arguments
parser = argparse.ArgumentParser(description="Dock a protein-ligand pair with AutoDock vina")
parser.add_argument("-p", help="Protein identifier (e.g. 6AZ1)", default="6RXC")
parser.add_argument("-l", help="Ligand or drug identifier (e.g. DB11989)", default="DB00336")
parser.add_argument("-v", help="Path to AutoDock Vina executable", default="~/Downloads/autodock_vina_1_1_2_linux_x86/bin/vina")
args = parser.parse_args()

if not os.path.exists("./Docking-Results"):
    os.mkdir("./Docking-Results")

VINA_PATH = args.v
PDBQT_DIR = "./PDBQTs"
CONFIG_DIR = "./CONFIGs"
LIGANDS_DIR = "./LIGANDs"
OUTPUT_DIR = "./Docking-Results"

# Get ligands and proteins
protein = args.p
ligand = args.l

receptor_file = f"{PDBQT_DIR}/{protein}sc.clean.pdbqt"
ligand_file = f"{LIGANDS_DIR}/{ligand}.pdbqt"

# Get Vina configuration file
vina_config_file = f"{CONFIG_DIR}/{protein}sc_vina.config"

# Get Vina output files
log_file = f"{OUTPUT_DIR}/{protein}_{ligand}_vina_results.log"
out_file = f"{OUTPUT_DIR}/{protein}_{ligand}_vina_results.pdbqt"
command = f"{VINA_PATH} --receptor {receptor_file} --ligand {ligand_file} --config {vina_config_file} --out {out_file} --log {log_file}"
os.system(command)

# Add top pose score to Docking Score Report
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
    os.system(f'echo "{protein},{ligand},{score}" >> Docking-Score-Report.csv')
    
except:
    pass

print("FINISHED WITH SUCCESS!")