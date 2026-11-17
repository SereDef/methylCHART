import os
import sys
import json
import warnings
import contextlib
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from pcntoolkit import (NormData, 
                        BLR, # HBR
                        BsplineBasisFunction, LinearBasisFunction, 
                        NormativeModel) #, plot_centiles, plot_qq

# https://www.ejwagenmakers.com/
# https://github.com/opherdonchin 
# https://www.biorxiv.org/content/10.64898/2026.02.17.706268v2.abstract

seed = 73
proj_dir = "/home/s.defina/methylCHART"

# ============================================================================
# --- config ---
# ============================================================================

# Settings grid: one SLURM array index == one combination
# ============================================================================
BASE = {
  'algorithm':   'BLR',         # 'BLR' | 'HBR' (only BLR implemented currently)
  'split':       (0.7, 0.3), 
  'sample':      'multicohort', # 'GenR' | 'ALSPAC' | 'multicohort'
  'knot_method': None,   # None = quantile for Age/Age2, uniform for logAge2
  'cell_comp':   False,  # include UniLIFE fractions (linear, k-1)
  'likelihood':  None,   # None = SHASH for beta, Gauss for Mval; or 'SHASH' | 'Gauss'
  } 

COMBOS = [
    {**BASE, 
    'data_scale': s, 
    'age_transform': a }
    for s in ('beta', 'Mval')
    for a in ('Age2', 'logAge2')
]

idx = int(sys.argv[1]) if len(sys.argv) > 1 else int(os.environ.get('SLURM_ARRAY_TASK_ID', 0))
if idx >= len(COMBOS):
    sys.exit(f"idx {idx}: only {len(COMBOS)} combos -> use --array=0-{len(COMBOS) - 1}")

model_desc = COMBOS[idx]

# model_desc = {
#     'algorithm':     'BLR',        # 'BLR' | 'HBR'
#     'split':         (0.7, 0.3),
#     'sample':        'multicohort',# 'GenR' | 'ALSPAC' | 'multicohort'
#     'age_transform': 'Age2',       # 'Age' | 'Age2' | 'logAge2'
#     'data_scale':    'Mval',       # 'beta' | 'Mval'
#     'likelihood':    'Gauss',      # 'SHASH' | 'Gauss' (SHASH~beta, Gauss~Mval)
#     'cell_comp':      True,        # include UniLIFE fractions (linear, k-1)
# }

age_raw = model_desc['age_transform']
knot_method = model_desc['knot_method'] or ('quantile' if age_raw in ('Age', 'Age2') else 'uniform')
likelihood  = model_desc['likelihood']  or ('SHASH' if model_desc['data_scale'] == 'beta' else 'Gauss')


model_name = "_".join([
  model_desc['algorithm'], 
  likelihood,
  model_desc['data_scale'],
  age_raw, f"{knot_method}Knots",
  'uniLIFE' if model_desc['cell_comp'] else 'nocell',
  model_desc['sample']
])

save_dir = f"{proj_dir}/results/{model_name}"
log_dir = f"{save_dir}/logs"
os.makedirs(log_dir, exist_ok=True)
print(f"[{idx}] {model_name}\n[{idx}] {model_desc}")

# ============================================================================
# Read data & assemble model
# ============================================================================

# Read input data
dnam  = pd.read_csv(f"{proj_dir}/data/dnam_data.csv", index_col=0)
pheno = pd.read_csv(f"{proj_dir}/data/pheno_data.csv", index_col=0)

data = pd.merge(dnam, pheno, on='Sample_ID') # Should be already in the same order 

# ---- Select sample -------------------------------------------------------

if model_desc['sample'] in ('GenR', 'ALSPAC'):
    key = {'GenR': 'GENR', 'ALSPAC': 'ALSPAC'}[model_desc['sample']]
    data = data.loc[data['Cohort'] == key]
  
print(pd.crosstab(data['Array'], data['Cohort']))

# ---- Select response -------------------------------------------------------

# Should not be any NAs - handled in previous stages
cpg_with_na = dnam.columns[dnam.isna().any()].tolist()
response_vars = sorted(set(dnam.columns) - set(['Sample_ID'] + cpg_with_na))

if model_desc['data_scale'] == 'Mval':
    eps = 1e-3
    b = np.clip(data[response_vars].to_numpy(), eps, 1 - eps)
    data[response_vars] = np.log2(b / (1 - b))

# ---- Select covariates -------------------------------------------------------

# Age form 
if age_raw == 'logAge2':
    assert (data['Age2'] > -1).all(), "log1p needs age > -1 (check preterm Age2 coding)"
    data['logAge2'] = np.log1p(data['Age2'])
    age_col = 'logAge2'
else:
    age_col = age_raw
    
# Cell types 
cell_cols = [c for c in pheno.columns if c.startswith('cell_')] # 19 UniLIFE fractions

# Keep age first and from one (first) cell
covariates = [age_col] + (cell_cols[:-1] if model_desc['cell_comp'] else []) 

# ---- Select batches  -------------------------------------------------------

batch_effects = ["Sex", "Array"] + (
  ["Cohort"] if model_desc['sample'] == 'multicohort' else []) # Period? + Plate

# ============================================================================
# Create NormData & split
# ============================================================================

norm_data = NormData.from_dataframe(
    name="cpgdata", dataframe=data, subject_ids='Sample_ID', # IDC?
    covariates=covariates, 
    batch_effects=batch_effects,
    response_vars=response_vars, 
    remove_outliers=False, #?
    z_threshold=10,
    # remove_Nan=True,
)

# Inspect
# norm_data.coords
# norm_data.data_vars

train, test = norm_data.train_test_split(
    splits=model_desc['split'], split_names=["train", "test"], random_state=seed)

# Inspect
df_train = train.to_dataframe()
df_test  = test.to_dataframe()


# ============================================================================
# Define model template
# ============================================================================

Algo = BLR  # if model_desc['algorithm'] == 'BLR' else HBR

model_template = Algo(
    name = model_name,
    
    # --- Covariates ---
    # B-spline on mean age (column 0); cell fractions pass through linearly
    basis_function_mean = BsplineBasisFunction(basis_column=0, 
                                               nknots=5,knot_method=knot_method,
                                               degree=2),
    # basis_function_var=LinearBasisFunction(),  # linear by default..?
    # model variance / noise is a function of the age + cell-types
    heteroskedastic = True,
    
    # --- Likelihood ---
    # other opthion 'WarpLog', 'WarpBoxCox', 'WarpAffine', 'WarpCompose'
    warp_name = ('WarpSinhArcsinh' if likelihood == 'SHASH' else None),
    # warp_reparam=True, # use a reparameterized warp function ??
    
    # --- Batch model ---
    # i.e. (Sex, Array, Cohort) 
    fixed_effect=True,  # model offsets in the mean for each batch
    fixed_effect_slope=True,        # model a fixed effect in the slope of the mean for each batch
    fixed_effect_slope_indices=[0], # Batch-specific slope on age only, not on cells
    fixed_effect_var=False,       # model a fixed effect in the intercept of the variance
    fixed_effect_var_slope=False, # model a fixed effect in the slope of the variance for each batch
    
    # n_iter=300,
    # tol=1e-4,
    # ard=False, # automatic relevance determination
    # l_bfgs_b_l = 0.1,
)

def model_config(template, model_name, base_dir=proj_dir):
  return NormativeModel(template_regression_model=template,
      savemodel=False, # for sharing and trasfering or to avoid refitting
      evaluate_model=True, # model fit metrics
      saveresults=True, # per-subject Z logp and centiles
      saveplots=True,
      save_dir=save_dir,
      inscaler="standardize", # "minmax", "robminmax", or "none"
      outscaler="standardize")

model = model_config(template=model_template, model_name=model_name)


# ============================================================================
# Fit CpGs in parallel
# ============================================================================
N_WORKERS = int(os.environ.get("SLURM_CPUS_PER_TASK", 10))
N_CHUNKS = min(N_WORKERS, len(response_vars))  # a wider stride would give empty chunks

def fit_chunk(i):
    """Fit CpGs i, i+N, i+2N, ... Forked worker, so train/test/model are inherited
    copy-on-write and nothing is pickled. Striding spreads slow CpGs evenly.
    Workers share save_dir: pcntoolkit merges the result CSVs under a FileLock."""
    tr = train.isel(response_vars=slice(i, None, N_CHUNKS))
    te = test.isel(response_vars=slice(i, None, N_CHUNKS))
    with open(f"{log_dir}/chunk_{i:03d}.log", "w", buffering=1) as log, \
            contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        warnings.simplefilter("always")  # every convergence warning, not just the first
        model.fit_predict(tr, te)
    return i, int(tr.response_vars.size)


if __name__ == "__main__":
    print(f"[{idx}] {len(response_vars)} CpGs | {N_WORKERS} workers | {N_CHUNKS} chunks")
    failed = []
    with ProcessPoolExecutor(max_workers=N_WORKERS, mp_context=mp.get_context("fork")) as ex:
        futures = {ex.submit(fit_chunk, i): i for i in range(N_CHUNKS)}
        for fut in as_completed(futures):
            try:
                i, n = fut.result()
                print(f"  chunk {i:03d} ok ({n} CpGs)", flush=True)
            except Exception as e:
                failed.append(futures[fut])
                print(f"  chunk {futures[fut]:03d} FAILED: {type(e).__name__}: {e}", flush=True)

    # stitch the per-chunk logs into one file, keeping attribution
    with open(f"{save_dir}/run_log.txt", "w", encoding="utf-8") as out:
        for name in sorted(os.listdir(log_dir)):
            if name.startswith("chunk_") and name.endswith(".log"):
                out.write(f"\n===== {name} =====\n")
                with open(f"{log_dir}/{name}", encoding="utf-8") as f:
                    out.write(f.read())
    
    # the one shared file every worker writes, unprotected by a lock: check it survived
    if model.savemodel:  # one shared file workers write unprotected: check it survived
        with open(f"{save_dir}/model/normative_model.json", encoding="utf-8") as f:
            json.load(f)

    print(f"[{idx}] done -> {save_dir}/results/statistics_cpgdata.csv")
    if failed:
        sys.exit(f"[{idx}] FAILED chunks {failed} - see {log_dir}")

# modelfit = model.fit_predict(train, test)
# 
# # Show the evaluation metrics from the train / test set
# print(round(train.get_statistics_df().T, 5))
# print(round(test.get_statistics_df().T, 5))

# === More plotting ============================================================
# import os
# import matplotlib.pyplot as plt
# from pcntoolkit import plot_centiles_advanced
# 
# out_dir = f"{proj_dir}/results/{model_name}/more_plots"
# os.makedirs(out_dir, exist_ok=True)
# 
# if model_desc['sample'] == 'multicohort':
#   color_batch = {"Cohort": ["GENR"], 'Sex': ['girl','boy'], 'Array': ['450k','EPICv1','EPICv2']}
# else:
#   color_batch = {'Sex': ['girl','boy'], 'Array': ['450k']}
# 
# figs = plot_centiles_advanced(
#     model,
#     scatter_data=train,
#     centiles=[0.05, 0.5, 0.95],
#     # Highlight a sepcific batch
#     batch_effects=color_batch,
#     # Show other data not belonging to the groups above
#     show_other_data=True,
#     # Remove the batch effects from the data (what the data would have looked like if all data was from the same batch)
#     harmonize_data=False,
# )
# 
# for i, fig in enumerate(figs):
#     # Use the CpG name from the plot title if there is one, otherwise a number
#     title = fig.axes[0].get_title() if fig.axes else ""
#     cpg = next((w for w in title.replace(",", " ").split() if w.startswith("cg")), f"plot_{i:02d}")
#     fig.savefig(f"{out_dir}/{cpg}_highlight.png", dpi=150, bbox_inches="tight")
#     plt.close(fig)   # free memory, avoids the "More than 20 figures" warning
# 
