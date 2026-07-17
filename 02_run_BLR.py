import pandas as pd
from pcntoolkit import NormData, BLR, BsplineBasisFunction, NormativeModel, plot_centiles, plot_qq

# Read input data
data = pd.read_csv("~/methylCHART/data/DNAm_data.csv")

# Specify modes structure
covariates = ["Age"]
batch_effects = ["sex", "Array", "IDC"] # "Period"celltype unilife 
response_vars = ["cpg"]

# Create a NormData object from the dataframe
norm_data = NormData.from_dataframe(
    name="cpgdata",
    dataframe=data,
    covariates=covariates,
    batch_effects=batch_effects,
    response_vars=response_vars,
    remove_outliers=False,
    z_threshold=10,
)

# Inspect
norm_data.coords
norm_data.data_vars

# Split
train, test = norm_data.train_test_split() # default: 80/20

# Inspect
df_train = train.to_dataframe()
df_test  = test.to_dataframe()

# ================ BLR ======================

# Specify the regression model
basic_blr = BLR(
    name="bspline",
    # We use a B-spline basis expansion for the mean, so the predicted mean is a 
    # smooth function of the covariates
    basis_function_mean=BsplineBasisFunction(degree=3, nknots=5),
    # The variance is a function of the covariates
    heteroskedastic=True
)

batch_blr = BLR(
    name="batched",
    # We use a B-spline basis expansion for the mean
    basis_function_mean=BsplineBasisFunction(degree=3, nknots=5),
    # The variance is a function of the covariates
    heteroskedastic=True,
    # Model the batch effects (Sex, Array)
    fixed_effect=True, # Model offsets in the mean for each individual batch effect
    fixed_effect_slope=True, # Model fixed effect in the slope of the mean for each individual batch effect
    fixed_effect_var_slope=True, # Model fixed effect in the slope of the variance for each individual batch effect
)

warped_blr = BLR(
    name="warped",
    # We use a B-spline basis expansion for the mean
    basis_function_mean=BsplineBasisFunction(degree=3, nknots=5),
    # The variance is a function of the covariates
    heteroskedastic=True,
    # Configure a sinh-arcsinh warp
    warp_name="warpsinharcsinh", 
    # Model the batch effects (Sex, Array, Period)
    fixed_effect=True, # Model offsets in the mean for each individual batch effect
    fixed_effect_slope=True, # Model fixed effect in the slope of the mean for each individual batch effect
    fixed_effect_var_slope=True, # Model fixed effect in the slope of the variance for each individual batch effect
)

# Configure the normative model
def model_config(template, model_name, base_dir="/home/s.defina/methylCHART"):
  model = NormativeModel(
      template_regression_model=template, # we select our BLR model
      savemodel=False, # for sharing and trasfering or to avoid refitting
      evaluate_model=True, # model fit metrics
      saveresults=True, # per-subject Z logp and centiles
      saveplots=True,
      save_dir=f"{base_dir}/{model_name}",
      inscaler="standardize", # "minmax", "robminmax", or "none"
      outscaler="standardize")
  return model

model0 = model_config(template=basic_blr, model_name="BLR_simple")
model0fit = model0.fit_predict(train, test)

model1 = model_config(template=batch_blr, model_name="BLR_batches")
model1fit =model1.fit_predict(train, test)

model2 = model_config(template=warped_blr, model_name="BLR_warped")
model2fit = model2.fit_predict(train, test)


# Show the evaluation metrics from the train / test set
print(round(model2fit.get_statistics_df().T, 5))
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
