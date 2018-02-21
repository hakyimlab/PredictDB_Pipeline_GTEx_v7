# PredictDB_Pipeline_GTEx_v7
PredictDB Pipeline for GTEx v7 data

Comparison of V6p and V7 models can be found [here](http://hakyimlab.org/post/2017/v7-v6p-analysis/)

# On GTEX V8 usage:

Data was gathered from GTEx' dbGap exchange folder in 'provisional files'.
Genotype files were created from GTEx shapeit2 vcfs.
Expression files were taken from eqtl folder and sliced to use our individuals of interest. A similar thing was done for covariate files.
Snp annotation was sliced and diced from the dbSnp golden path.

For the scripts to work, the following folders must be created in `model_training`:
```bash
mkdir summary
mkdir weights
mkdir dbs    
mkdir analysis
mkdir covariances
```

Then, run in this order:

```
make_dbs.R
filter_dbs.R
create_covariances.py
```

These scripts need to be modified as they contain hardcoded paths and file names and patterns.