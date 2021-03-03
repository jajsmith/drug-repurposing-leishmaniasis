import pickle
import pandas as pd
import numpy as np
# from pubchempy import get_compounds
from tqdm import tqdm
import sys
from time import sleep

# ins = ['1EVZ', '1P33', '6RXC']
# pubchem_id， smile, score for each protein
# ins = [1, 2, 4, 10, 12, 19, 20, 21, 26, 27, 28, 29, 30, 31]
#
#
# mapping = pd.read_csv('data/cid_drug_smile.csv', index_col=0).to_numpy()
#
# # scoring = scoring[:100]
# for i in ins:
#     with open('out/tmp/save_folder_6AZ1_{}/output_list.pkl'.format(i),
#               'rb') as f:
#         scoring = pickle.load(f)
#     with tqdm(total=len(scoring), file=sys.stdout) as pbar:
#         f = open('6AZ1_{}_results.csv'.format(i), 'w')
#         f.write(','.join(['score', 'cid', 'smile']) + '\n')
#         for c, _, s in scoring:
#             c = int(c.split('_')[1])
#             f.write(','.join([s, mapping[c, 0], mapping[c, 1]]) + '\n')
#             pbar.set_description('processed-{}: '.format(i))
#             pbar.update(1)
#         f.close()

# merge results on proteins for each drug dataset
drugset = ['in-trial', 'drugcentral', 'endogenous', 'world']
protset = ['Q6TUJ5', 'O97193', 'O15826', 'Q4QHJ8', 'Q4QBL1', 'P48499', 'Q27686']

for i in [1]:
    out = []
    for j in range(7):
        with open(
                'MONN/{}/save_folder_{}/results_aggregation/output_list.pkl'.format(
                        drugset[i], protset[j]), 'rb') as f:
            out.append(pickle.load(f))
    out = pd.DataFrame(np.array(out).reshape((-1, 3)),
                       columns=['drug', 'prot', 'score'])
    out.sort_values(by=['score'], inplace=True)
    out['score'] = out['score'].astype(np.float)
    out.to_csv('MONN/DeepPurpose_{}_rank.csv'.format(drugset[i]), index=False)
    a = out['score'][out['drug'] == 'ZINC000022059930'][out['prot'] == 'Q6TUJ5']

    # compute correlation
    monn = pd.read_csv('MONN/preds_{}_ordered.csv'.format(drugset[i]),
                       usecols=[0, 1, 3])
    monn['Unnamed: 0'] = ['id-{}'.format(i) for i in monn['Unnamed: 0']]
    f = open('MONN/{}/merged_rank.csv'.format(drugset[i]), 'w')
    f.write("zinc_id,pid,Kd,DeepPurpose_score\n")
    for z, p, s in monn.to_numpy():
        tmp = out['score'][out['drug'] == z][out['prot'] == p].to_numpy()
        tmp = tmp[0] if tmp.tolist() else np.NAN
        f.write("{},{},{},{}\n".format(z, p, s, tmp))
    f.close()
