
## Goal

Build a reference model for DNA methylation across development

## Reference data 

* **Generation R**: `7,143` observations from `2533` individuals (5 time points)
* **ALSPAC**: `6,106` observations from `3155` individuals (8 time points)

Total: **`13,249`** observations from `5688` individuals

<ins>Pre-processing</ins>: functional normalization with/+ batch correction

<ins>*N* features</ins>: **`366,764`** (overlapping 450K, EPICv1 & EPICv2)

Input size: **`36.2 GB`** (16.7 ALSPAC + 19.5 GENR)

## Algorithm

* **BLR** with **sin-arcsinh (SHASH) likelihood** [empirical Bayes]
* (alternative) *HBR* with *SHASH* or *Beta* likelihood [full Bayesian estimation]


## Model specification

Covariates:

* **Age** 
  - with mean centered gestational age (`Age2`) or age at birth fixed at 0 (`Age`)
  - *Mean model*: B-spline expansion (5 knots, 3rd degree poly)
  - *Variance model*: heteroskedastic (linear)
  
* *Unilife cell-type composition* [TODO] (linear, k-1)

Batches: 

* **Cohort**
* **Array**
* **Sex**

```
fixed_effect=True,  # model offsets in the mean for each batch
fixed_effect_slope=True,      # model a fixed effect in the slope of the mean for each batch
fixed_effect_var=False,       # model a fixed effect in the intercept of the variance
fixed_effect_var_slope=False, # model a fixed effect in the slope of the variance for each batch
```

> Note: ignore repeated measures 

## Other choices / info

* *No outlier removal* (`z_threshold`= 10 - very lenient) 
* *70-30 split* (automatically stratified by batch)
* Benchmark: 28 CpGs (GenR only) **~3 min** and 81M of space
* env: `Python 3.12`, `pcntoolkit 1.3.0` 

* Benchmark HBR ? [TODO]
* Extract convergence warnings? 
  - (e.g. All-ill conditioned matrix, not pd, overflow encountered in dot, divide by zero encountered in log)
* Model comparison
* Smart chunking
* Optimize output? 