# %%
import enum
from os import replace
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors as ds
import numpy as np
import umap
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go

# %%
# select descriptor modules
modules = ['ExactMolWt', 'FpDensityMorgan1', 'FpDensityMorgan2', 'FpDensityMorgan3', 'HeavyAtomMolWt', 'MaxAbsPartialCharge', 'MaxPartialCharge', 'MinAbsPartialCharge', 'MinPartialCharge', 'MolWt', 'NumValenceElectrons']
# modules = ['MolWt', 'BertzCT', 'NumAliphaticRings', 'NumAromaticRings', 'NumHAcceptors', 'NumHDonors', 'NumHeteroatoms', 'NumRotatableBonds', 'MolLogP', 'fr_SH', 'RingCount', 'NumAromaticCarbocycles', 'NumAromaticHeterocycles', 'NOCount']

def generate_descriptor_features(smile_list, module_list) -> np.array:
    """
    Parameters:
        - smile_list/module_list: array-like formate
    """
    features = np.zeros((len(smile_list), len(modules)))
    for i, smile in enumerate(smile_list):
        mol = Chem.MolFromSmiles(smile)
        features[i] = [ds.__dict__[f](mol) for f in modules]

    return features
    
def generate_embedding(scaled_data, dim=3, mod='pca', seed = 20):
    if mod == 'umap':
        reducer = umap.UMAP(n_components=dim, random_state=seed)
        embeddings = reducer.fit_transform(scaled_data) 
    if mod == 'pca':
        pca = PCA(n_components=3, random_state=seed)
        pca.fit(scaled_data)
        embeddings = pca.transform(scaled_data)
    
    return embeddings, mod

# %%
#********************************************************
# visualization of screening of in-trial drugs with 5RXC
#********************************************************

scores = pd.read_csv('scores_with_smiles.tsv', sep='\t', usecols=[0, 2, 3])

# generate features
features = generate_descriptor_features(scores['SMILES'], modules)

# dimensional reduction of features
random_seed = 18
df = pd.DataFrame(features, columns=modules)
df['zincID'] = scores['In-Trial Drug ZINC ID']
df['Indaba_score'] = scores['INDABA']
df = df.dropna()

data = df[modules].values
scaled_data = StandardScaler().fit_transform(df[modules].values)

# construct data for plot
embeddings, mod = generate_embedding(mod='pca', seed=random_seed)
embed_data = df[['zincID', 'Indaba_score']]
embed_data['x'] = embeddings[:, 0]
embed_data['y'] = embeddings[:, 1]
embed_data['z'] = embeddings[:, 2]

# plot
fig = px.scatter_3d(embed_data, x='x', y='y', z='z', color='Indaba_score', color_continuous_scale='ice', hover_name='zincID', title=f'with random-seed: {random_seed}')
# fig.write_html(f'{mod}_inTrial_6RXC_seed_{random_seed}.html')

# %%
#*******************************************************
# find those with the smallest euclidean distance in 3D to the top drugs
#*******************************************************

map_idx_to_smile = {}
# process drugCentral
drug_central = pd.read_csv('drugCentral.csv', usecols=[1, 4])
for smile, idx in drug_central.values:
    map_idx_to_smile[f'drugCentral_{idx}'] = smile
# process intrils
intrial = pd.read_csv('in-trials.csv')
for idx, smile in intrial.values:
    map_idx_to_smile[idx] = smile

#%%
indices = list(map_idx_to_smile.keys())
smiles = [map_idx_to_smile[i] for i in indices]

# generate features
features = generate_descriptor_features(smiles, modules)
features = StandardScaler().fit_transform(features)



#%%
# read top pairs
top_score_25 = pd.read_csv('top_pairs_25.tsv', sep='\t')
top_drugs = top_score_25['molecule'].unique()

indices = np.array(indices)
top_features = np.zeros((len(top_drugs), len(modules)))
for i, mol in enumerate(top_drugs):
    top_features[i] = features[np.where(indices==mol)][0]

#%%
# compute distance matrix
dist_matrix = np.sqrt(np.sum((top_features[None, :] - features[:, None])**2, -1))
dist = pd.DataFrame(dist_matrix, columns=top_drugs, index=indices)

#%%
df = pd.DataFrame(features, columns=modules)
df['molecule'] = indices
df = df.dropna()

#%%
embeddings, mod = generate_embedding(df[modules], mod='umap', seed=8)
top_embeddings = np.zeros((len(top_drugs), 3))
for i, mol in enumerate(top_drugs):
    top_embeddings[i] = embeddings[np.where(df['molecule'] == mol)][0]

fig = go.Figure()

fig.add_trace(go.Scatter3d(
    x=embeddings[:, 0],
    y=embeddings[:, 1],
    z=embeddings[:, 2],
    mode='markers',
    marker=dict(size=2),
    hovertext=df['molecule'],
    name='in-trial + drug central'
))
fig.add_trace(
    go.Scatter3d(
        x=top_embeddings[:, 0],
        y=top_embeddings[:, 1],
        z=top_embeddings[:, 2],
        mode='markers',
        marker=dict(color=2, line_width=1),
        hovertext=top_drugs,
        name='top drug'
    )
)
fig.write_html('tmp.html')


# %%
