import pandas as pd
from openbabel import pybel
import os

## Get in-trials drugs
print("Getting in-trial drug files")

df = pd.read_csv('in-trials.csv')
identifiers = list(df['zinc_id'])
smiles = list(df['smiles'])

if not os.path.exists("../Ligands"):
    os.mkdir("../Ligands")

for i in range(len(identifiers)):
    smile = smiles[i]
    identifier = identifiers[i]
    try:
        mol = pybel.readstring("smi", smile)
        mol.addh()
        mol.localopt(steps=500)
        mol.write("pdbqt", "../Ligands/{}.pdbqt".format(identifier), overwrite=True)
        mol.write("sdf", "../Ligands/{}.sdf".format(identifier), overwrite=True)
    except:
        continue

## Get DrugCentral drugs

print("Getting Drug Central drug files")

df = pd.read_csv('drugCentral.csv', index_col=0)
identifiers = ["DC" + str(i) for i in list(df['ID'])]
smiles = list(df['SMILES'])

for i in range(len(identifiers)):
    smile = smiles[i]
    identifier = identifiers[i]
    try:
        mol = pybel.readstring("smi", smile)
        mol.addh()
        mol.localopt(steps=500)
        mol.write("pdbqt", "../Ligands/{}.pdbqt".format(identifier), overwrite=True)
        mol.write("sdf", "../Ligands/{}.sdf".format(identifier), overwrite=True)
    except:
        continue

## Get DrugBank drugs

print("Getting DrugBank drug files")

df = pd.read_csv('drugBank_leishmania.csv')
identifiers = list(df['DRUGBANK_ID'])
smiles = list(df['SMILES'])

for i in range(len(identifiers)):
    smile = smiles[i]
    identifier = identifiers[i]
    try:
        mol = pybel.readstring("smi", smile)
        mol.addh()
        mol.localopt(steps=500)
        mol.write("pdbqt", "../Ligands/{}.pdbqt".format(identifier), overwrite=True)
        mol.write("sdf", "../Ligands/{}.sdf".format(identifier), overwrite=True)
    except:
        continue

print("Finished with success")