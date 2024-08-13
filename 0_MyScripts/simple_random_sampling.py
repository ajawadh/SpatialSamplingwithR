import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import numpy as np
import os

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

print(pd_grdVoorst_df.head())

n = 40
N = len(pd_grdVoorst_df)

np.random.seed(314)
units = np.random.choice(N, size=n, replace=True)

my_sample = pd_grdVoorst_df.iloc[units, :]

print(my_sample.head())

cellsize = 25


def jitter(values, amount):
    return values + np.random.uniform(-amount, amount)


my_sample.loc[:, "s1"] = jitter(my_sample["s1"], amount=cellsize / 2)
my_sample.loc[:, "s2"] = jitter(my_sample["s2"], amount=cellsize / 2)

print(my_sample.head())
