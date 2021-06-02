# DeepPurpose Protein-Drug Binding Prediction

For running protein-drug binding predictions on Leishmania targets.

## Requirements

Tested on Ubuntu 18.04

Load the conda environment file:

```
cd DeepPurpose/
conda env create -f environment.yml
conda activate deeppurpose
```

## Generating Predictions

### Step 1. Running Predictions

```
python run.py
```

### Step 2. Final Output Ranking

```
python j.py
```