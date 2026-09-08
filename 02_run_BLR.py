import pandas as pd
from pcntoolkit import NormData, BLR, BsplineBasisFunction, LinearBasisFunction, NormativeModel #, plot_centiles, plot_qq

seed = 73

# Read input data
cpg_data = pd.read_csv(f"~/methylCHART/data/DNAm_data.csv", index_col=0)
cov_data = pd.read_csv(f"~/methylCHART/data/cov_data.csv", index_col=0)

data = pd.merge(cov_data, cpg_data, on='Sample_ID')

cpg_with_na = cpg_data.columns[cpg_data.isna().any()].tolist()
other_cpgs = list(set(cpg_data.columns) - set(['Sample_ID']+cpg_with_na))

# Specify modes structure
covariates = ["Age"]
batch_effects = ["sex", "Array"] # "IDC" "Period"celltype unilife 
response_vars = other_cpgs

# do not account for for idc, rather downstream in the z scores -- velocity and what not
# split 50-50 - but why do I need to use the test set with external data? 
# birth age what to do? Bing remind to put in contact with N
# https://www.ejwagenmakers.com/
# https://github.com/opherdonchin 
# https://www.biorxiv.org/content/10.64898/2026.02.17.706268v2.abstract

# Create a NormData object from the dataframe
norm_data = NormData.from_dataframe(
    name="cpgdata",
    dataframe=data,
    subject_ids='Sample_ID',
    covariates=covariates,
    batch_effects=batch_effects,
    response_vars=response_vars,
    remove_outliers=False,
    z_threshold=10,
    
    # remove_Nan=True,
)

# Inspect
norm_data.coords
norm_data.data_vars

# Split
train, test = norm_data.train_test_split(0.5, random_state=seed)

# Inspect
df_train = train.to_dataframe()
df_test  = test.to_dataframe()

# ================ BLR ======================

warped_blr = BLR(
    name="warpedBLR",
    # We use a B-spline basis expansion for the mean
    basis_function_mean=BsplineBasisFunction(degree=3, nknots=5),
    # basis_function_var=LinearBasisFunction(),  # linear by default?
    # The variance is a function of the covariates
    heteroskedastic=True, # noise modeling
    # Configure a sinh-arcsinh warp
    warp_name="WarpSinhArcsinh", 
    warp_reparam=True, # use a reparameterized warp function
    # Model the batch effects (Sex, Array, Period)
    fixed_effect=True, # Model offsets in the mean for each batch effect
    fixed_effect_var=True, 
    fixed_effect_slope=False, # Model fixed effect in the slope of the mean for each batch effect
    fixed_effect_var_slope=False, # Model fixed effect in the slope of the variance for each batch effect
    n_iter=300,
    tol=1e-4,
    ard=False, # automatic relevance determination
    l_bfgs_b_l = 0.1,
)

# Configure the normative model
def model_config(template, model_name, base_dir="/home/s.defina/methylCHART"):
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

model = model_config(template=warped_blr, model_name="warpedBLR")
modelfit = model.fit_predict(train, test)


# Show the evaluation metrics from the train / test set
print(round(train.get_statistics_df().T, 5))
print(round(test.get_statistics_df().T, 5))


# Cast back to a dataframe 
# new_df = train.sel(response_vars="cpg").to_dataframe()
# new_df.head()

# figs = plot_centiles(model, scatter_data=train)
# qcfs = plot_qq(test, plot_id_line=True)
# 
# # Save the figure as PNG
# def save_plots(plots, name, loc = "/home/s.defina/methylCHART"):
#   for i, p in enumerate(plots, start = 1):
#       p.savefig(f"{loc}/{name}{i}.png", dpi = 300, bbox_inches = "tight")
# 
# save_plots(figs, 'plot_blr_centile')
# save_plots(qcfs, 'plot_blr_qq')

# Create a BLR model with heteroskedastic noise
# model = NormativeModel(BLR(heteroskedastic=True),
#                     inscaler='standardize',
#                     outscaler='standardize')
# 
# model.fit_predict(train, test)
