import pickle
import numpy as np
import matplotlib.pyplot as plt

with open('out/output_list_0.pkl', 'rb') as f:
    result_0 = np.array(pickle.load(f))

with open('out/output_list_1.pkl', 'rb') as f:
    result_1 = np.array(pickle.load(f))

print()


# ######################## PROTEIN_0 ########################
score_0 = result_0[:, 2].astype(np.float)
print(score_0.mean())                           # mean: 3,858.82
print(score_0.min())                            # min:  0.06
print(score_0.max())                            # max:  3,673,380.03

# ######################## PROTEIN_1 ########################
score_1 = result_1[:, 2].astype(np.float)
print(score_1.mean())                           # mean: 44,163.14
print(score_1.min())                            # min:  0.98
print(score_1.max())                            # max:  5,543,577.84


plt.figure()
n = 10
plt.hist(score_0[score_0 < n], bins=1000, label='Prot_0')
plt.hist(-1 * score_1[score_1 < n], bins=1000, label='Prot_1')
plt.legend()
plt.show()


# ######################## PROTEIN_ALL ########################
n = 100
result_all = np.concatenate((result_0, result_1), axis=0)
sorted_result_all = result_all[result_all[:, 2].astype(np.float).argsort()]
topn = sorted_result_all[:n, :]



import pandas as pd
chem_map = pd.read_csv("data/chem_all.csv").to_numpy()
topn_idx = topn[:, 0].astype(np.str)
topn_idx_1 = result_1[:n, 0].astype(np.str)

topn_smile = [chem_map[int(idx.split('_')[1]), 0] for idx in topn_idx]
top20_smile_1 = [chem_map[int(idx.split('_')[1]), 0] for idx in topn_idx_1]






from rdkit import Chem
from rdkit.DataStructs import FingerprintSimilarity

sm = np.zeros((n, n))
for i in range(n):
    for j in range(i, n):
        m1, m2 = Chem.MolFromSmiles(top20_smile_1[i]), Chem.MolFromSmiles(top20_smile_1[j])
        sm[i, j] = FingerprintSimilarity(Chem.RDKFingerprint(m1), Chem.RDKFingerprint(m2))

sm = sm + sm.T - np.eye(n)

from sklearn.cluster import AffinityPropagation

af = AffinityPropagation().fit(sm)
cluster_centers_indices = af.cluster_centers_indices_
labels = af.labels_

n_clusters_ = len(cluster_centers_indices)
af.cluster_centers_

print('Estimated number of clusters: %d' % n_clusters_)
# print("Silhouette Coefficient: %0.3f"
#       % metrics.silhouette_score(X, labels, metric='sqeuclidean'))

sent_0 = ['CHEM_6372', 'CHEM_92102', 'CHEM_16957', 'CHEM_31681', 'CHEM_92515',
          'CHEM_87689', 'CHEM_40322', 'CHEM_16498', 'CHEM_4109', 'CHEM_16518']

out = {}
out['overlapping_chem'] = ['CHEM_10758','CHEM_2724','CHEM_3777','CHEM_62242','CHEM_6372','CHEM_73001','CHEM_89869']
out['protein_0'] = {'AAs': 'MTAPTVPVALVTGAAKRLGRSIAEGLHAEGYAVCLHYHRSAAEANALSATLNARRPNSAITVQADLSNVATAPVSGADGSAPVTLFTRCAELVAACYTHWGRCDVLVNNASSFYPTPLLRNDEDGHEPCVGDREAMETATADLFGSNAIAPYFLIKAFAHRFAGTPAKHRGTNYSIINMVDAMTNQPLLGYTIYTMAKGALEGLTRSAALELAPLQIRVNGVGPGLSVLVDDMPPAVWEGHRSKVPLYQRDSSAAEVSDVVIFLCSSKAKYITGTCVKVDGGYSLTRA',
                    'picked_chem': sent_0,
                    'top100_cluster_label': [0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,
                                             1,  1,  1,  1,  1,  1,  1,  1,  1,  9,  2,  2,  2,  2,  2,  2,  2,
                                             2,  1,  3,  9,  3,  3,  3,  3,  3,  3,  5,  9,  9,  0,  5,  9,  9,
        9,  9,  9,  5,  4,  4,  4,  4,  5,  5,  5,  6,  6,  6,  6,  7,  7,
        7,  7,  7,  7,  5,  8,  8,  8,  8,  8,  8,  8,  8,  5,  9,  9,  9,
        5,  9,  9,  5,  5,  9, 10, 10, 10, 10, 10, 10, 10, 10,  1],
                    'cluster_centers_indices': [ 0, 20, 27, 39, 55, 59, 63, 66, 73, 82, 92]}
out['protein_1'] = {'AAs': 'MSRAAARFKIPMPETKADFAFPSLRAFSIVVALDMQHGIGDGESIPWRVPEDMTFFKNQTTLLRNKKPPTEKKRNAVVMGRKTWESVPVKFRPLKGRLNIVLSSKATVEELLAPLPEGQRAAAAQDVVVVNGGLAEALRLLARPLYCSSIETAYCVGGAQVYADAMLSPCIEKLQEVYLTRIYATAPACTRFFPFPPENAATAWDLASSQGRRKSEAEGLEFEICKYVPRNHEERQYLELIDRIMKTGIVKEDRTGVGTISLFGAQMRFSLRDNRLPLLTTKRVFWRGVCEELLWFLRGETSAQLLADKDIHIWDGNGSREFLDSRGLTENKEMDLGPVYGFQWRHFGADYKGFEANYDGEGVDQIKLIVETIKTNPNDRRLLVTAWNPCALQKMALPPCHLLAQFYVNTDTSELSCMLYQRSCDMGLGVPFNIASYALLTILIAKATGLRPGELVHTLGDAHVYRNHVDALKAQLERVPHAFPTLIFKEERQYLEDYELTDMEVIDYVPHPAIKMEMAV',
                    'picked_chem': 'To Be decided',
                    'top100_cluster_label': [ 0, 10, 10,  1, 10,  2,  2,  2,  2,  2,  2,  2, 10,  1,  1,  1,  1,
        9,  1,  1,  1,  1,  1, 10,  3,  3,  3,  3,  4,  4,  4,  4,  4,  5,
        5,  5,  5,  5,  5,  5,  5,  4,  4,  4,  4,  4,  4,  0,  0,  0,  1,
        1,  1,  1,  1, 10,  6,  6,  6,  6, 10, 10, 10,  8, 10, 10,  0, 10,
       10, 10,  7,  7,  7,  7,  7,  7,  7,  7,  8,  8,  8,  9,  9,  9,  4,
        0, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11],
                    'cluster_centers_indices': [ 0,  3,  8, 27, 32, 35, 57, 71, 78, 82, 91, 99]}
out['protein_0']['rank_aff_score'] = result_0.tolist()
out['protein_1']['rank_aff_score'] = result_1.tolist()

mapping = {v:k for k, v in chem_map}
out['chem_to_smail'] = mapping

import pickle

with open('info.pkl', 'wb') as f:
    pickle.dump(out, f)
