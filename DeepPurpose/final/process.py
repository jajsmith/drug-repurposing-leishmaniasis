import numpy as np
import pandas as pd

google = pd.read_csv("final/google_sequence.csv").drop_duplicates()
rcsb = pd.read_csv('final/RCSB_sequence.csv').drop_duplicates()

mapping = {}

seq_ids = []
for seq in google['seq']:
    if mapping.get(seq) is None:
        mapping[seq] = len(mapping)
    seq_ids.append(mapping[seq])
google['seq_id'] = [f"seq-{i}" for i in seq_ids]
google = google[['pdb', 'seq_id', 'seq']]

seq_ids = []
for seq in rcsb['seq']:
    if mapping.get(seq) is None:
        mapping[seq] = len(mapping)
    seq_ids.append(mapping[seq])
rcsb['seq_id'] = [f"seq-{i}" for i in seq_ids]
rcsb = rcsb[['pdb', 'seq_id', 'seq']]

rcsb[['seq_id', 'seq']].to_csv('final/r_seq.csv', index=False)
google[['seq_id', 'seq']].to_csv('final/g_seq.csv', index=False)

rcsb.to_csv("rcsb_pdb-seqID-seq.csv", index=False)
google.to_csv("google_pdb-seqID-seq.csv", index=False)

# ##################################
# prepare data for virtual screening
prot = pd.read_csv('final/g_seq.csv').to_numpy()
drug = pd.read_csv('final/chem.csv').to_numpy()

out = []
# f = open('rcsb.tsv', 'w')
# f.write(f'SmileID\tSmile\tTargetID\tTarget Sequence\n')
for seqID, seq in prot:
    for smileID, smile in drug:
        out.append([smileID, smile, seqID, seq])
        # f.write(f'{smileID}\t{smile}\t{seqID}\t{seq}\n')
# f.close()
out = pd.DataFrame(out, columns=['SmileID', 'Smile', 'TargetID', 'TargetSequence'])
out.to_csv('google.tsv', sep='\t', index=False)


# ####################################
import pickle
with open('google/results_aggregation/output_list_VS.pkl', 'rb') as f:
    result = pickle.load(f)
data = pd.DataFrame(result, columns=['smileID', 'seqID', 'score'])
data.to_csv('google.csv', index=False)


