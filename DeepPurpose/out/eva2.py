import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.DataStructs import FingerprintSimilarity

#           path                                           binding score
results = ["out/tmp/save_folder_1EVZ_1/output_list.pkl",        # 0.22
           "out/tmp/save_folder_1P33_1/output_list.pkl",        # 0.18
           "out/tmp/save_folder_1T10_1/output_list.pkl",        # 1.08
           "out/tmp/save_folder_3HQQ_1/output_list.pkl",        # 0.71
           "out/tmp/save_folder_1AMK_1/output_list.pkl"]        # 2.65

# #################################################
result_path = results[1]
# #################################################

with open(result_path, 'rb') as f:
    result = np.array(pickle.load(f))

score = result[:, 2].astype(np.float)
print("mean: {}".format(score.mean()))
print("min: {}".format(score.min()))
print("max: {}".format(score.max()))

plt.figure()
n = 100
plt.hist(score[:n], bins=10)
plt.title("Hist for top {} --".format(n) + result_path.split('/')[2].split('_')[2])
plt.legend()
plt.show()


topn = result[:n, :]
topn_idx = topn[:, 0].astype(np.str)

chem_map = pd.read_csv("data/chem_all.csv").to_numpy()
topn_smile = [chem_map[int(idx.split('_')[1]), 0] for idx in topn_idx]


sm = np.zeros((n, n))
for i in range(n):
    for j in range(i, n):
        m1, m2 = Chem.MolFromSmiles(topn_smile[i]), Chem.MolFromSmiles(topn_smile[j])
        sm[i, j] = FingerprintSimilarity(Chem.RDKFingerprint(m1), Chem.RDKFingerprint(m2))

sm = sm + sm.T - np.eye(n)

from sklearn.cluster import AffinityPropagation

af = AffinityPropagation().fit(sm)
cluster_centers_indices = af.cluster_centers_indices_
labels = af.labels_

n_clusters_ = len(cluster_centers_indices)
print("{} clusters: ".format(n_clusters_))
print("    Center index: {}".format(cluster_centers_indices.tolist()))
print("    Labels: {}".format(labels.tolist()))

send = {'1EVZ': ['CHEM_833', 'CHEM_84524', 'CHEM_6372', 'CHEM_28096', 'CHEM_16023'],
        '1P33': ['CHEM_6372', 'CHEM_40322', 'CHEM_4109', 'CHEM_8472'],
        '3HQQ': ['CHEM_6372', 'CHEM_40322', 'CHEM_4109', 'CHEM_16498'],
        '1T10': ['CHEM_6372', 'CHEM_40322', 'CHEM_3777', 'CHEM_38064', 'CHEM_74497']}


