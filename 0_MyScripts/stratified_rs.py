import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
import numpy as np
import os
import seaborn as sns
import scipy.stats as stats
import matplotlib
import pandas as pd

matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import statsmodels.stats.api as sms
import statsmodels.api as sm

# import pandas as pd

# os.chdir()
print(os.getcwd)

# Enable automatic conversion of R objects to pandas DataFrames
pandas2ri.activate()

# Load the .rda file
robjects.r["load"]("data/grdVoorst.rda")

# Assuming the R data frame is called 'grdVoorst'
r_grdVoorst_df = robjects.r["grdVoorst"]

# Convert to pandas DataFrame if it's an R DataFrame
if isinstance(r_grdVoorst_df, robjects.vectors.DataFrame):
    with localconverter(robjects.default_converter + pandas2ri.converter):
        pd_grdVoorst_df = pandas2ri.rpy2py(r_grdVoorst_df)

    print(type(pd_grdVoorst_df))  # Should output <class 'pandas.core.frame.DataFrame'>
    print(pd_grdVoorst_df.columns)  # Now you can access the columns
else:
    print(f"The object is not a DataFrame, it's a {type(df)}")

print(pd_grdVoorst_df.head())
# Count the number of grids for each stratum
N_h = pd_grdVoorst_df["stratum"].value_counts()
print(N_h.sum())

# Calculate the proportion of each stratum
w_h = N_h / N_h.sum()

# Define the total sample size
n = 40

# Calculate the number of samples for each stratum
n_h = np.round(n * w_h).astype(int)

# Print the result
print(n_h)

n_h[n_h == n_h.max()] = n_h[n_h == n_h.max()] - 1
print(n_h)


# Define a function to perform stratified sampling with replacement
def stratified_sample(
    df, strata_col, stratum_sample_size, strata_weight, prop_strata_unit, sizes
):
    df_stratified = pd.DataFrame()

    for stratum, size in sizes.items():
        print(stratum)
        print(size)
        stratum_df = df[df[strata_col] == stratum]
        stratum_df = stratum_df.copy()
        prob = 1 - float((1 - (1 / len(stratum_df)))) ** float(size)
        stratum_df.loc[:, strata_weight] = stratum_df.shape[0] / df.shape[0]
        stratum_df.loc[:, prop_strata_unit] = prob
        stratum_df.loc[:, stratum_sample_size] = size
        sample = stratum_df.sample(n=size, replace=True, random_state=314)
        df_stratified = pd.concat([df_stratified, sample], axis=0)

    return df_stratified


# Sort the strata to ensure consistent ordering
ord = N_h.index.sort_values()
print(ord)
print(n_h[ord].to_dict())

# Perform the stratified sampling
mysample = stratified_sample(
    pd_grdVoorst_df,
    "stratum",
    "stratum_sample_size",
    "stratum_weight",
    "prop_strata_unit",
    n_h[ord].to_dict(),
)

# Apply jitter to columns 's1' and 's2'
# mysample.loc[:, "s1"]  = mysample["s1"] + np.random.uniform(-12.5, 12.5, size=len(mysample))
# mysample.loc[:, "s2"] = mysample["s2"] + np.random.uniform(-12.5, 12.5, size=len(mysample))


def jitter(values, amount):
    return values + np.random.uniform(-amount, amount, size=values.shape)


cellsize = 25
mysample.loc[:, "s1"] = jitter(mysample["s1"], amount=cellsize / 2)
mysample.loc[:, "s2"] = jitter(mysample["s2"], amount=cellsize / 2)

# View the sampled data
print(mysample)

# Calculate the estimator of the mean of `z` within each stratum
mz_h = mysample.groupby("stratum")["z"].mean()

# Calculate the overall estimated weighted mean
mz = (w_h * mz_h).sum()

# Calculate the estimated variance of `z` within each stratum
S2z_h = mysample.groupby("stratum")["z"].var(ddof=1)

# Calculate the estimated sampling variance of the estimator of stratum means
v_mz_h = S2z_h / n_h

# Calculate the estimated standard error of the estimator weighted mean
se_mz = np.sqrt((w_h**2 * v_mz_h).sum())

# Print the results
print(f"Weighted Mean (mz): {mz}")
print(f"Standard Error (se_mz): {se_mz}")

# Compile all estimates into a dictionary
estimates = {
    "stratum_estimates": {
        stratum: {
            "mean": mz_h[stratum],
            "variance": S2z_h[stratum],
            "sampling_variance": v_mz_h[stratum],
            "n": n_h[stratum],
        }
        for stratum in mz_h.index
    },
    "population_estimates": {
        "weighted_mean": mz,
        "weighted_standard_error": se_mz,
        "n": n,
    },
}

# Print the compiled estimates dictionary
print(estimates)

# Convert stratum estimates to a DataFrame
stratum_df = pd.DataFrame.from_dict(estimates["stratum_estimates"], orient="index")

# Convert population estimates to a DataFrame
population_df = pd.DataFrame(estimates["population_estimates"], index=["Population"])

# Combine both DataFrames
combined_df = pd.concat([stratum_df, population_df])

# Display the final DataFrame
print(combined_df)

# estimating the weighted population mean directly using the inclusion probabilities "prob_strata_unit"
mz_estimator = np.sum(mysample["z"] / mysample["prop_strata_unit"]) / N_h.sum()
print(mz_estimator)


# Calculate the unique probability for each stratum
pi_h = mysample.groupby("stratum")["prop_strata_unit"].unique().apply(lambda x: x[0])
print(pi_h)
print(pi_h * N_h)

# Compute the sum of the product of unique probabilities and total units in each stratum
result = (pi_h * N_h).sum()

# Print the result, we see that the reverse engineered sample size from inclusion probability estimates is not exactly 40 due to rounding errors
print(result)

# Estimating the mean with simple unweighted mean
mz_unweighted = np.mean(mysample["z"])
print(mz_unweighted)
# this unweighted mean is slightly biased. The bias would have been substantially
#  larger if an equal number of units would have been selected from each stratum,
# leading to much larger differences in the inclusion probabilities among the strata.
