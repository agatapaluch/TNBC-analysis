import pandas as pd

df_id = pd.read_csv("COSMIC_v3.4_ID_GRCh37.txt", sep='\t')
df_dbs = pd.read_csv("COSMIC_v3.4_DBS_GRCh38_exome.txt", sep='\t')

# IMMUCAN
# ID: Signatures to keep (10% / >=5 samples): ['ID6', 'ID8', 'ID9']
# DBS: Signatures to keep (10% / >=5 samples): ['DBS17']

keep_columns_immucan_id = ['Type', 'ID6', 'ID8', 'ID9']
keep_columns_immucan_dbs = ['Type', 'DBS17']


# SYNERGY
# ID: Signatures to keep (10% / >=10 samples): ['ID6', 'ID8', 'ID9']
# DBS: Signatures to keep (10% / >=10 samples): ['DBS2', 'DBS17']

keep_columns_synergy_id = ['Type', 'ID6', 'ID8', 'ID9']
keep_columns_synergy_dbs = ['Type','DBS2', 'DBS17']

df_immucan_id = df_id.loc[:, keep_columns_immucan_id]
df_immucan_dbs = df_dbs.loc[:, keep_columns_immucan_dbs]

df_synergy_id = df_id.loc[:, keep_columns_synergy_id]
df_synergy_dbs = df_dbs.loc[:, keep_columns_synergy_dbs]

df_immucan_id.to_csv('IMMUcan_CUSTOM_COSMIC_v3.4_ID_GRCh37.txt', sep='\t', index=False)
df_immucan_dbs.to_csv('IMMUcan_CUSTOM_COSMIC_v3.4_DBS_GRCh38_exome.txt', sep='\t', index=False)

df_synergy_id.to_csv('SYNERGY_CUSTOM_COSMIC_v3.4_ID_GRCh37.txt', sep='\t', index=False)
df_synergy_dbs.to_csv('SYNERGY_CUSTOM_COSMIC_v3.4_DBS_GRCh38_exome.txt', sep='\t', index=False)


