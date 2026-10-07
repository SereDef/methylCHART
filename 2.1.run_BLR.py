import pandas as pd
from pcntoolkit import NormData, BLR, BsplineBasisFunction, LinearBasisFunction, NormativeModel #, plot_centiles, plot_qq

# https://www.ejwagenmakers.com/
# https://github.com/opherdonchin 
# https://www.biorxiv.org/content/10.64898/2026.02.17.706268v2.abstract

seed = 73

proj_dir = "/home/s.defina/methylCHART"

model_desc = {
  'algorithm': 'BLR', # HBR
  'split': (0.7, 0.3), # or other
  'likelihood': 'SHASH',
  'sample': 'multicohort',
  'birth_age_model': 'Age' # Mean centered gestational age or Age for birth fixed at 0
}

model_name = f"{model_desc['algorithm']}_{model_desc['likelihood']}_{model_desc['sample']}_{model_desc['birth_age_model']}"

# Read input data
dnam  = pd.read_csv(f"{proj_dir}/data/dnam_data.csv", index_col=0)
pheno = pd.read_csv(f"{proj_dir}/data/pheno_data.csv", index_col=0)

data = pd.merge(dnam, pheno, on='Sample_ID') # Should be already in the same order 

# Select sample 
if model_desc['sample'] == 'GenR':
  data = data.loc[data['Cohort'] == 'GENR']
  
elif model_desc['sample'] == 'ALSPAC':
  data = data.loc[data['Cohort'] == 'ALSPAC']
  
print(pd.crosstab(data['Array'], data['Cohort']))

# Should not be any NAs - handled in previous stages
cpg_with_na = dnam.columns[dnam.isna().any()].tolist()
complete_obs_cpgs = list(set(dnam.columns) - set(['Sample_ID']+cpg_with_na))

# ===== Specify model structure ================================================
covariates = [model_desc['birth_age_model']] # Age or Age2

batch_effects = ["Sex", "Array", *(
  ["Cohort"] if model_desc['sample']=='multicohort' else [])] # "IDC" "Period"celltype unilife 

response_vars = complete_obs_cpgs

# === Preprocess dataset =======================================================
norm_data = NormData.from_dataframe(
    name="cpgdata",
    dataframe=data,
    subject_ids='Sample_ID', # IDC?
    covariates=covariates,
    batch_effects=batch_effects,
    response_vars=response_vars,
    remove_outliers=False, # ?
    z_threshold=10,
    # remove_Nan=True,
)

# Inspect
# norm_data.coords
# norm_data.data_vars

# Split
train, test = norm_data.train_test_split(splits=model_desc['split'], 
    split_names=["train", "test"], random_state=seed)

# Inspect
df_train = train.to_dataframe()
df_test  = test.to_dataframe()

# === Define model template ====================================================


model_template = (BLR if model_desc['algorithm']=='BLR' else HBR)(
    name = model_name,
    # --- Age model ---
    # use a B-spline basis expansion for the mean
    basis_function_mean = BsplineBasisFunction(degree=3, nknots=5),
    # basis_function_var=LinearBasisFunction(),  # linear by default..?
    # model variance / noise is a function of the age
    heteroskedastic = True,
    # --- Likelihood ---
    # other opthion 'WarpLog', 'WarpBoxCox', 'WarpAffine', 'WarpCompose'
    warp_name = 'WarpSinhArcsinh' if model_desc['likelihood']=='SHASH' else None, 
    # warp_reparam=True, # use a reparameterized warp function ??
    # --- Batch model ---
    # i.e. (Sex, Array, Cohort) 
    fixed_effect=True,  # model offsets in the mean for each batch
    fixed_effect_slope=True,      # model a fixed effect in the slope of the mean for each batch
    fixed_effect_var=False,       # model a fixed effect in the intercept of the variance
    fixed_effect_var_slope=False, # model a fixed effect in the slope of the variance for each batch
    # n_iter=300,
    # tol=1e-4,
    # ard=False, # automatic relevance determination
    # l_bfgs_b_l = 0.1,
)

# Configure the normative model
def model_config(template, model_name, base_dir=proj_dir):
  model = NormativeModel(
      template_regression_model=template, # we select our BLR model
      savemodel=False, # for sharing and trasfering or to avoid refitting
      evaluate_model=True, # model fit metrics
      saveresults=True, # per-subject Z logp and centiles
      saveplots=True,
      save_dir=f"{base_dir}/results/{model_name}",
      inscaler="standardize", # "minmax", "robminmax", or "none"
      outscaler="standardize")
  return model

model = model_config(template=model_template, model_name=model_name)
modelfit = model.fit_predict(train, test)

# Show the evaluation metrics from the train / test set
print(round(train.get_statistics_df().T, 5))
print(round(test.get_statistics_df().T, 5))

# === More plotting ============================================================
import os
import matplotlib.pyplot as plt
from pcntoolkit import plot_centiles_advanced

out_dir = f"{proj_dir}/results/{model_name}/more_plots"
os.makedirs(out_dir, exist_ok=True)

if model_desc['sample'] == 'multicohort':
  color_batch = {"Cohort": ["GENR"], 'Sex': ['girl','boy'], 'Array': ['450k','EPICv1','EPICv2']}
else:
  color_batch = {'Sex': ['girl','boy'], 'Array': ['450k']}

figs = plot_centiles_advanced(
    model,
    scatter_data=train,
    centiles=[0.05, 0.5, 0.95],
    # Highlight a sepcific batch
    batch_effects=color_batch,
    # Show other data not belonging to the groups above
    show_other_data=True,
    # Remove the batch effects from the data (what the data would have looked like if all data was from the same batch)
    harmonize_data=False,
)

for i, fig in enumerate(figs):
    # Use the CpG name from the plot title if there is one, otherwise a number
    title = fig.axes[0].get_title() if fig.axes else ""
    cpg = next((w for w in title.replace(",", " ").split() if w.startswith("cg")), f"plot_{i:02d}")
    fig.savefig(f"{out_dir}/{cpg}_highlight.png", dpi=150, bbox_inches="tight")
    plt.close(fig)   # free memory, avoids the "More than 20 figures" warning

