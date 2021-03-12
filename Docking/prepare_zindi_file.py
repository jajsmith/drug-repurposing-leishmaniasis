"""
Arguments

1 input directory
2 input.pdbqt file
3 output folder
4 protein PDB file
"""
import os
import sys

input_dir = sys.argv[1]
input_file = sys.argv[2]
output_dir = sys.argv[3]
protein_pdb = sys.argv[4]

zindi_file = input_file.replace(".pdbqt", ".pdb")

# Step 1 - Convert Vina docking PDBQT result to raw PDB
os.system(f"obabel -ipdbqt -f 1 -l 1 {input_dir}/{input_file} -opdb > {output_dir}/tmp.pdb")

# Step 2 - Clean up "ATOM" lines
os.system(f"sed 's/ATOM  /HETATM/g' {output_dir}/tmp.pdb > {output_dir}/cleaned_tmp.pdb")
os.system(f'grep "HETATM" {output_dir}/cleaned_tmp.pdb > {output_dir}/cleaned_ligand.pdb')
os.system(f'grep "CONECT" {output_dir}/cleaned_tmp.pdb >> {output_dir}/cleaned_ligand.pdb')
os.system(f'grep "MASTER" {output_dir}/cleaned_tmp.pdb >> {output_dir}/cleaned_ligand.pdb')
os.system(f'grep "END" {output_dir}/cleaned_tmp.pdb >> {output_dir}/cleaned_ligand.pdb')
# Step 3 - Concatenate protein and ligand file
os.system(f"cat {protein_pdb} {output_dir}/cleaned_ligand.pdb > {output_dir}/{zindi_file}")
