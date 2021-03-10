import requests
import os
import urllib.request
import pandas as pd
from Bio import SeqIO

df = pd.read_csv('leishmania_pdbs_in_RCSB.csv')
pdb_ids = list(df['PDB ID'])
pdb_ids = sorted(list(set(pdb_ids)))

print("Downloading original PDB files from RCSB...")

if not os.path.exists("./PDBs"):
    os.mkdir("PDBs")

count = 0
for pdb_id in pdb_ids:
    try:
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb_id}.pdb", f"PDBs/{pdb_id}.pdb")
        count += 1
    except:
        continue

print(f"A total of {count} PDB files were downloaded")

#%% Write single chain PDB files

print("Generating single chain PDB files...")

def lineBelongsToChain(line, chain="A"):
    words = line.split()
    line_type = words[0]
    if line_type == "ATOM":
        return words[4] == chain
    elif line_type == "TER":
        return words[3] == chain
    else:
        return False

original_pdbs = os.listdir("PDBs")
original_pdbs = [file for file in original_pdbs if ".pdb" in file]

if not os.path.exists("./single-chain-PDBs"):
    os.mkdir("single-chain-PDBs")

for file in original_pdbs:
    f = open("PDBs/" + file, 'r')
    new_file = "single-chain-PDBs/" + file.split('.')[0] + "sc.clean.pdb"
    g = open(new_file, 'w')
    line = f.readline()
    while len(line) > 0:
        if lineBelongsToChain(line):
            g.write(line)
        line = f.readline()
        continue
    f.close()
    g.close()

#%% Extract sequences from each PDB file downloaded

print("Extracting amino acid sequences from PDB files...")

def getSequenceFromPDB(filename):
    handle = SeqIO.parse(filename, "pdb-atom")
    record = list(handle)[0]
    seq = str(record.seq)
    return seq

protein_files = os.listdir("single-chain-PDBs")
print(f"There are {len(protein_files)} single chain protein PDB files!")

if not os.path.exists("./Sequences"):
    os.mkdir("Sequences")

f = open("Sequences/pdb_sequences.fa", 'w')

print("Writing FASTA file with all sequences...")

for file in protein_files:
    filename = f"single-chain-PDBs/{file}"
    try:
        seq = getSequenceFromPDB(filename)
    except:
        continue
    header = f">{file}"
    f.write(header + '\n' + seq + '\n')

f.close()

#%% Generate PDBQT files of each protein

print("Generating PDBQT files...")

if not os.path.exists("./PDBQTs"):
    os.mkdir("PDBQTs")

protein_files = os.listdir("single-chain-PDBs")

for file in protein_files:
    filename = f"single-chain-PDBs/{file}"
    
    # Step 1 - Get raw PDBQT file
    command_1 = f"obabel -ipdb {filename} -opdbqt > ./PDBQTs/tmp.pdbqt"
    os.system(command_1)

    # Step 2 - Clean up protein's PDBQT file
    pdbqt_file = "PDBQTs/" + file.replace('.pdb','.pdbqt')
    command_2 = f'grep "ATOM" ./PDBQTs/tmp.pdbqt > {pdbqt_file}'
    os.system(command_2)

#%% Create Configuration files for AutoDock vina

print("Creating AutoDock Vina configuration files")

if not os.path.exists("./CONFIGs"):
    os.mkdir("CONFIGs")

# Function to retrieve atom coordinates
def getSearchBoxCoordinatesFromPDB(pdb_file):
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure('PROT', pdb_file)
    model = structure[0]
    chains = list(model.child_dict.keys())
    chain = model[chains[0]]
    atoms = list(chain.get_atoms())
    coord_matrix = np.zeros([len(atoms), 3])
    for i,atom in enumerate(atoms):
        coord_matrix[i] = atom.get_coord()
    results = {}
    results['mean_xyz'] = list(coord_matrix.mean(axis=0))
    results['min_xyz'] = list(coord_matrix.min(axis=0))
    results['max_xyz'] = list(coord_matrix.max(axis=0))
    return results

# Function to retrieve Vina search box parameters
def getSearchBoxConfigParams(search_box_coordinates):
    min_x, min_y, min_z = search_box_coordinates['min_xyz']
    max_x, max_y, max_z = search_box_coordinates['max_xyz']
    dist_x = abs(max_x - min_x) # distance between min and max
    dist_y = abs(max_y - min_y)
    dist_z = abs(max_z - min_z)
    center_x = min_x + (dist_x/2)
    center_y = min_y + (dist_y/2)
    center_z = min_z + (dist_z/2)
    size_x = dist_x * 1.25
    size_y = dist_y * 1.25
    size_z = dist_z * 1.25
    results = {}
    results['center_x'] = round(center_x,2)
    results['center_y'] = round(center_y,2)
    results['center_z'] = round(center_z,2)
    results['size_x'] = round(size_x)
    results['size_y'] = round(size_y)
    results['size_z'] = round(size_z)
    return results

# Function to write a Vina config file given box parameters
def writeVinaConfigFile(filename, box_params, modes=10, cpu=5):
    center_x = box_params['center_x']
    center_y = box_params['center_y']
    center_z = box_params['center_z']
    size_x = box_params['size_x']
    size_y = box_params['size_y']
    size_z = box_params['size_z']
    template = f"""center_x = {center_x}
center_y = {center_y}
center_z = {center_z}

size_x = {size_x}
size_y = {size_y}
size_z = {size_z}

num_modes = {modes}
cpu = {cpu}
seed= 8

exhaustiveness = 8"""

    f = open(filename, 'w')
    f.write(template)
    f.close()
    print(f"Finished writing {filename}!")

# Read arguments
input_pdb_folder = "single-chain-PDBs"
output_config_folder = "CONFIGs"

pdb_files = os.listdir(input_pdb_folder)
pdb_files = [file for file in pdb_files if '.pdb' in file]

for input_pdb_file in pdb_files:
    try:
        protein_name = input_pdb_file[0:input_pdb_file.index('.')]

        # Get Box parameters
        search_box_coordinates = getSearchBoxCoordinatesFromPDB(f"{input_pdb_folder}/{input_pdb_file}")
        box_params = getSearchBoxConfigParams(search_box_coordinates)

        # Write config file
        output_file = f"{output_config_folder}/{protein_name}_vina.config"
        writeVinaConfigFile(output_file, box_params)
    except:
        continue

#%% Finish run
print("Finished with success")