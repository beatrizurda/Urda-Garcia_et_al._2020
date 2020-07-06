# Urda-Garcia_et_al._2020
Code for the manuscript 'From comorbidities to gene expression fingerprints and back'

## Data
Poner un [enlace](https://pip.pypa.io/en/stable/) a GEO y GREIN.

## Steps and Code
### Analyzing disease similarities based on gene expression at the disease level
1. Run the RNA-seq pipeline at the disease level
```bash
run_rnaseq_pipeline_for_disease.R
```
2. Perform functional enrichment and produce it's heatmap-based visualization
```bash
molecular_insight_heatmap.R
```
5. Generate the networks
- Obtain the distance between all pairwise comparisons
```bash
build_disease_level_network.py
```
- Generate the networks using the predefined distances
```bash
generating_networks.R
```
6. Compute the overlap of our networks with the health electronic records-derived network from Hidalgo et al.
```bash
generating_networks.R
```
7. Run the entire anlysis and compute the overlap at the icd-9 level
- Merge the initial data to group diseases that correspond to the same ICD9 code `generating_SE_objects_icd9_level.R`
- Run the RNA-seq pipeline for each ICD9 code `prepare_pipeline_runs.R` & `run_rnaseq_pipeline_for_disease.R`
- Compute the distances between each ICD9 pairs `build_ICD_level_network.R`
- Generate the networks `generating_networks.R`
- Compute the overlap `network_overlap_icd.R`


### Analyzing disease similarities based on gene expression at the meta-patient level
1. Generate the meta-patients
```bash
defining_meta_patients.R
```
2. Run the RNA-seq pipeline at the meta-patient level
```bash
DEanalysis_for_metapatients.R
```
3. Build the Stratified Similarity Network (SSN) that connects meta-patients and diseases
- Obtain the distance between all pairwise comparisons
```bash
build_metapatient_disease_network.py 
```
- Generate the networks using the predefined distances
```bash
generating_networks.R
```

### Differential Variability Analysis
1. Perform Differential Variability Analysis for each disease
```bash
exploring_variance_distribution_disease_level.R 
```
2. Analyze the correlation of the metric Distance to Median (DM) with sample size and average gene expression.
```bash
DM_behaviour.R
```
3. Perform functional enrichment and produce it's heatmap-based visualization
```bash
molecular_insight_heatmap.R
```

