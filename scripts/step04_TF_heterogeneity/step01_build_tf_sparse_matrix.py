# Replace <LOCAL_DATA_ROOT>/<SERVER_DATA_ROOT> with your own paths.
# Build sparse matrix
import pandas as pd

# Read test.txt (assumed to be tab-delimited)
df = pd.read_csv("<LOCAL_DATA_ROOT>/36-AMD-retry/promoter-retry/ref_alt_TF.txt", sep="\t")

# Get Ref and Alt sets
ref_set = set(df["Ref"].dropna())
alt_set = set(df["Alt"].dropna())

# Use union as row names
all_items = sorted(ref_set | alt_set)

# Build matrix
matrix_df = pd.DataFrame(0, index=all_items, columns=["Ref", "Alt"])
matrix_df.loc[list(ref_set), "Ref"] = 1
matrix_df.loc[list(alt_set), "Alt"] = 1

# Output to txt file (tab-delimited)
matrix_df.to_csv("<LOCAL_DATA_ROOT>/36-AMD-retry/promoter-retry/ref_alt_TF_sparse_matrix.txt", sep="\t")

print("Sparse matrix file generated")
print(matrix_df.head())
