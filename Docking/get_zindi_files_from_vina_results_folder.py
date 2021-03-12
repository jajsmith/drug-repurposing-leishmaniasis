import os

if not os.path.exists("./Zindi-Ready-Files"):
    os.mkdir("./Zindi-Ready-Files")

vina_dir = "./Docking-Results"
PDB_DIR = "./single-chain-PDBs"
output_dir = "./Zindi-Ready-Files"

# Get ligand names
vina_files = os.listdir(vina_dir)
vina_files = [file for file in vina_files if '.pdbqt' in file]
proteins = [file.split("_")[0] for file in vina_files]

for i in range(len(vina_files)):
    pdbqt_file = vina_files[i]
    pdb_file = f"{PDB_DIR}/{proteins[i]}sc.clean.pdb"
    print(f"Converting file {pdbqt_file}")
    command = f"python prepare_zindi_file.py {vina_dir} {pdbqt_file} {output_dir} {pdb_file}"
    os.system(command)

print("Finished with success!")
