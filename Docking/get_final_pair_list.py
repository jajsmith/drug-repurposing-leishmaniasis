import pickle
import pandas as pd
import numpy as np
import sys
from time import sleep
import os
from Bio import SeqIO

# DeepPurpose Scores

dp_scores = pd.read_csv('../DeepPurpose/final_dp.csv')

n = 5
good_scores = dp_scores.head(int(len(dp_scores)*(n/100)))
bad_scores = dp_scores.tail(int(len(dp_scores)*((100-n)/100)))


proteins = good_scores['prot_seq'].drop_duplicates()
drugs = good_scores['drug_smiles'].drop_duplicates()

# Drug Selection

df_a = pd.read_csv('./indaba_drugs/drugCentral.csv', index_col=0)
identifiers_a = ["DC" + str(i) for i in list(df_a['ID'])]
smiles_a = list(df_a['SMILES'])

df_b = pd.read_csv('./indaba_drugs/in-trials.csv')
identifiers_b = list(df_b['zinc_id'])
smiles_b = list(df_b['smiles'])

df_c = pd.read_csv('./indaba_drugs/drugBank_leishmania.csv')
identifiers_c = list(df_c['DRUGBANK_ID'])
smiles_c = list(df_c['SMILES'])

final_drugs = []
for i, d in enumerate(smiles_a):
    if d in drugs.values:
        final_drugs.append(identifiers_a[i])

for i, d in enumerate(smiles_b):
    if d in drugs.values:
        final_drugs.append(identifiers_b[i])

final_drugs = final_drugs + identifiers_c


# Protein Selection

input_file = './Sequences/pdb_sequences.fa'

fasta_sequences = SeqIO.parse(open(input_file),'fasta')

pdbs = []
seqs = []
for fasta in fasta_sequences:
    pdbs.append(fasta.id)
    seqs.append(fasta.seq)

final_proteins = []

for i, s in enumerate(seqs):
    if str(s) in proteins.values:
        pdbid = pdbs[i].strip('sc.clean.pdb')
        if not pdbid.startswith("A0"):
            final_proteins.append(pdbid)

            
# Final pair construction

pairs = []
for p in final_proteins:
    for d in final_drugs:
        pairs.append(pd.Series([p,d]))

df = pd.concat(pairs, axis=1).transpose()

print("Final Pairs: ", len(df))
print("Writinig to file...")

df.to_csv('final.csv', index=False, header=False)

print ("Done!")