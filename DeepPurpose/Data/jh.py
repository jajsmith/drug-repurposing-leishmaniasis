import pickle
import pandas as pd
import numpy as np
from pubchempy import get_compounds
from tqdm import tqdm
import sys
from time import sleep

ins = ['1EVZ', '1P33', '6RXC']
# pubchem_id， smile, score for each protein

chem_map = pd.read_csv('chem_all.csv', usecols=[0]).to_numpy().flatten().tolist()
# cid_map = []

# ccc = 0

start, end = int(sys.argv[1]), int(sys.argv[2])
with tqdm(total=end-start, file=sys.stdout) as pbar:
    f = open('{}-{}.csv'.format(str(start), str(end)), 'w')
    f.write(','.join(['DP_index', 'cid', 'smiles']) + '\n')
    for idx in range(start, end):
        c = chem_map[idx]
        cid = None
        try:
            cid = get_compounds(c, 'smiles')[0].cid
        except:
            print(idx)
        f.write(','.join([str(idx), str(cid), c]) + '\n')
        pbar.set_description('processed: %d' % (1 + idx))
        pbar.update(1)
        # sleep(0.2)
    f.close()



#
# with open('out/tmp/save_folder_{}_1/output_list.pkl'.format(ins[0]), 'rb') as f:
#     inn = pickle.load(f)
# inn = np.array(inn)
#
# smiles, pubchem_id, dp_score = [], [], []
# for c, _, s in inn:
#     c = chem_map[int(c.split('_')[-1])]
#     smiles.append(cid_map[c])
#






