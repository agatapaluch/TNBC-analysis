import pandas as pd

#df = pd.read_csv("DBS/DBS_IMMUcan_before_filtering/Assignment_Solution/Activities/Assignment_Solution_Activities.txt", sep='\t')
#df = pd.read_csv("DBS/DBS_SYNERGY_before_filtering/Assignment_Solution/Activities/Assignment_Solution_Activities.txt", sep='\t')

df = pd.read_csv("ID/ID_IMMUcan_before_filtering/Assignment_Solution/Activities/Assignment_Solution_Activities.txt", sep='\t')
df = pd.read_csv("ID/ID_SYNERGY_before_filtering/Assignment_Solution/Activities/Assignment_Solution_Activities.txt", sep='\t')


sample_col = "Samples"
sig_cols = [c for c in df.columns if c != sample_col]

# Sum of mutations in each row
row_sums = df[sig_cols].sum(axis=1)
# Exposure
exposure = df[sig_cols].div(row_sums, axis=0)

### THRESHOLD = 5%
threshold = 0.05

# In how many samples exposure > 5%
n_samples_5 = (exposure > threshold).sum(axis=0)
#print(n_samples_5)

# Keep only signatures which are active in >=5 samples
keep_5 = n_samples_5[n_samples_5 >= 5].index.tolist()
filtered_5 = df[[sample_col] + keep_5]

print("Signatures to keep (5% / >=5 samples):")
print(keep_5)

### THRESHOLD = 10%
threshold = 0.10

n_samples_10 = (exposure > threshold).sum(axis=0)
keep_10 = n_samples_10[n_samples_10 >= 10].index.tolist()
filtered_10 = df[[sample_col] + keep_10]

print("Signatures to keep (10% / >=10 samples):")
print(keep_10)

keep_5 = n_samples_10[n_samples_10 >= 5].index.tolist()
print("Signatures to keep (10% / >=5 samples):")
print(keep_5)

