import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import numpy as np
import os
import seaborn as sns
import matplotlib
import scipy.stats as stats

matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import statsmodels.stats.api as sms

# import pandas as pd

print(os.getcwd)

# Enable automatic conversion of R objects to pandas DataFrames
pandas2ri.activate()

# Load the .rda file
robjects.r["load"]("data/grdVoorst.rda")

# Assuming the R data frame is called 'df'
r_grdVoorst_df = robjects.r["grdVoorst"]

# Convert to pandas DataFrame if it's an R DataFrame
if isinstance(r_grdVoorst_df, robjects.vectors.DataFrame):
    with localconverter(robjects.default_converter + pandas2ri.converter):
        pd_grdVoorst_df = pandas2ri.rpy2py(r_grdVoorst_df)

    print(type(pd_grdVoorst_df))  # Should output <class 'pandas.core.frame.DataFrame'>
    print(pd_grdVoorst_df.columns)  # Now you can access the columns
else:
    print(f"The object is not a DataFrame, it's a {type(df)}")

# A dataframe with the centers of the selected cells of the discretisation grid,
# referred to as discretisation points.
# To void restricting the sampling units to the discretisation points,
# a simple random sample of points is selected in two stages. First,
# n times a grid cell is selected by SRSWR.
# Second, everytime a grid cell is selected, one point is selected fully
# randomly from that grid cell. This account for the infinite number of points
# in the population. This second step of this selection procedure is implemented
# with function "jitter". it adds random noise to the spatial coordinates of grid centers
# by drawing a continous uniform distribution unif(-c,c) with c half the side length of
# square grid cells. Compared to the above selection, here we respect that the
# population is actually infinite.

# SRSWOR from discretized grids with points at the center, each square is 25 m in length
n = 40
N = len(pd_grdVoorst_df)  # grdVoorst is organized in a discritized grid in a tibble
# where each square center is saved. this is the sampling frame.

np.random.seed(314)
units = np.random.choice(N, size=n, replace=True)
# python unit index generation is different and thus yields different results compared to R sampling.

""" R units = [
    6478,
    1408,
    1084,
    6531,
    5076,
    3491,
    3318,
    5830,
    332,
    1064,
    6362,
    1304,
    2525,
    5717,
    4705,
    3283,
    678,
    6282,
    5314,
    2276,
    4969,
    5483,
    6580,
    1175,
    6005,
    4214,
    7147,
    5535,
    2757,
    4636,
    5391,
    397,
    1874,
    4084,
    7437,
    3331,
    5065,
    6958,
    1732,
    3907,
]"""

# Check sampled indices
print("Sampled indices:", units)

my_sample = pd_grdVoorst_df.iloc[units, :]

print(my_sample.head())

cellsize = 25


def jitter(values, amount):
    return values + np.random.uniform(-amount, amount, size=values.shape)


my_sample.loc[:, "s1"] = jitter(my_sample["s1"], amount=cellsize / 2)
my_sample.loc[:, "s2"] = jitter(my_sample["s2"], amount=cellsize / 2)

print(my_sample.head())

# calculating the mean of sample
mz = np.mean(my_sample.z)

# multiplying this by the area alone is not useful (different units)
# calculated the volumen of soil in the grids for the first 0-30 cm
# each grid is multiplied by the number of grids N * size of grid * thikness of soil
# total SOC is then estimated by volume * bulk density of soil * SOC mean
# this is multiplied by 10-6 to obtain total mass of SOC in Mg (1000 kg)

vol_soil = N * 0.3 * 25**2
bd = 1500
tz = vol_soil * bd * mz * 10**-6
print(tz)

# ideally, one should also measure the bulk density at the sampling points
# thus, calculating the volumetric SOC concentration in kg m-3
# the estimated population mean of this volumetric SOM can then be multiplied by
# the total volume of soil in the study area

# if you sample 10000 times, the sampling distrbution of the horvitz-thompson
# estimator of the mean SOC. this estimator distribution of the mean is normally distributed
# compared to the actual population frequency distribution.
# if you would repeat the sampling infinite number of times
# and after scaling, the area under the curve would be 1
# and you would obtain the sampling distribution of the estimator
# of the population mean.

# a good estimator should be unbiased and with low variance.
# the Mean Squared Error (MSE) combines both: sum of variance and
# squared bias-

# Create an ECDF plot with a step function
fig, ax = plt.subplots(1, 2, figsize=(12, 6))

sns.ecdfplot(data=my_sample, x="z", stat="proportion", complementary=False, ax=ax[1])
sns.kdeplot(data=my_sample, x="z", ax=ax[0])
sns.histplot(data=my_sample, x="z", stat="density", binwidth=10, ax=ax[0])

# Customize the axes
plt.xlabel("SOM")
plt.ylabel("F")

# Display the plot
# plt.show()

# calculation of quantiles
# method 1:
z = my_sample["z"]


# Implementing Type 4 quantile calculation manually
def type4_quantile(data, probs):
    data = np.sort(data)
    n = len(data)
    quantiles = [(np.floor(n * p + 0.5) - 1) for p in probs]
    return [data[int(q)] for q in quantiles]


# Calculate quantiles at 0.25, 0.5, and 0.75
quantiles = type4_quantile(z, [0.25, 0.5, 0.75])

# Round the results to one decimal place
rounded_quantiles = np.round(quantiles, 1)

print(rounded_quantiles)

# method 2: using numpy
# Calculate quantiles at 0.25, 0.5, and 0.75
quantiles = np.quantile(z, [0.25, 0.5, 0.75])

# Round the results to one decimal place
rounded_quantiles = np.round(quantiles, 1)

print(rounded_quantiles)

# method 3: using pandas
# Calculate quantiles at 0.25, 0.5, and 0.75
quantiles = z.quantile([0.25, 0.5, 0.75])

# Round the results to one decimal place
rounded_quantiles = quantiles.round(1)

print(rounded_quantiles)


# computing the confidence interval of quantile estimates
def bootstrap_ci(data, q, alpha=0.05, n_bootstraps=10000):
    """
    Compute the confidence interval for a quantile using bootstrap resampling.

    Parameters:
    - data: 1D array-like, the data to calculate the quantile for.
    - q: float, the quantile to compute (e.g., 0.5 for median).
    - alpha: float, significance level (e.g., 0.05 for 95% confidence interval).
    - n_bootstraps: int, number of bootstrap samples.

    Returns:
    - Tuple (lower_bound, upper_bound) representing the confidence interval.
    """
    bootstrapped_quantiles = []

    for _ in range(n_bootstraps):
        sample = np.random.choice(data, size=len(data), replace=True)
        bootstrapped_quantiles.append(np.quantile(sample, q))

    lower_bound = np.percentile(bootstrapped_quantiles, 100 * (alpha / 2))
    upper_bound = np.percentile(bootstrapped_quantiles, 100 * (1 - alpha / 2))

    return lower_bound, upper_bound


# Example usage:
# Assuming mysample is a pandas DataFrame and 'z' is the column of interest
z = my_sample["z"]

# Compute the confidence interval for the 0.5 quantile (median) with alpha = 0.05
ci = bootstrap_ci(z, q=0.5, alpha=0.05)

print(f"Confidence interval for the median: {ci}")


def exact_binomial_ci(data, q, alpha=0.05):
    """
    Compute the exact binomial confidence interval for a quantile.

    Parameters:
    - data: 1D array-like, the data to calculate the quantile for.
    - q: float, the quantile to compute (e.g., 0.5 for median).
    - alpha: float, significance level (e.g., 0.05 for 95% confidence interval).

    Returns:
    - Tuple (lower_bound, upper_bound) representing the confidence interval.
    """
    n = len(data)
    k = int(np.floor(n * q))  # The rank of the quantile in the sorted data
    sorted_data = np.sort(data)

    # Exact binomial confidence limits
    ci_low = sms.binom_test(k, n, prop=q, alternative="two-sided") / 2
    ci_high = sms.binom_test(k + 1, n, prop=q, alternative="two-sided") / 2

    lower_bound_index = int(np.floor(n * ci_low))
    upper_bound_index = int(np.ceil(n * (1 - ci_high)))

    # Bound the indices within valid range
    lower_bound_index = max(0, lower_bound_index)
    upper_bound_index = min(n - 1, upper_bound_index)

    return sorted_data[lower_bound_index], sorted_data[upper_bound_index]


# Example usage:
# Assuming mysample is a pandas DataFrame and 'z' is the column of interest
z = my_sample["z"]

# Compute the confidence interval for the 0.5 quantile (median) with alpha = 0.05
ci = exact_binomial_ci(z, q=0.5, alpha=0.05)

print(f"Exact binomial confidence interval for the median: {ci}")
