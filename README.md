# Drug Repurposing: Leishmaniasis

Steps:

1. Load initial data into evaluation format
2. Run DeepPurpose on Protein-Compound pairs in `DeepPurpose/`
3. Rank Protein-Compound pairs by predicted binding affinity from DeepPurpose
4. Run Autodock Vina for top pairs in `Docking/`
5. Submit top predicted Protein-Compound pairs