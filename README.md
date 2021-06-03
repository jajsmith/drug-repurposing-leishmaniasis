# Drug Repurposing: Leishmaniasis

## Environment:

Use the provided conda environment in each directory when running code.

Tested on Ubuntu 18.04 with Conda.

## Steps:

1. Run the `Aggregation.ipynb` in `DeepPurpose/` to generate deepurpose binding scores for proteins and compounds. 
2. Run Autodock Vina for all pairs on top proteins in `Docking/`. Follow instructions in README to setup environment and run docking automatically for the `final.csv` set of Protein-Compound pairs.
3. Submit top predicted Protein-Compound pairs from `Docking/Zindi-Ready-Files`.