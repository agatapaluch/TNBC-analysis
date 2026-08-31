import pandas as pd

df = pd.read_csv("COSMIC_v3.4_SBS_GRCh38_exome.txt", sep='\t')

# IMMUCAN (non-ffpe repaired)
# Signatures to keep (10% / >=5 samples): ['SBS1', 'SBS2', 'SBS5', 'SBS13', 'SBS39', 'SBS87']

# SYNERGY (non-ffpe repaired)
#Signatures to keep (10% / >=5 samples): ['SBS1', 'SBS2', 'SBS3', 'SBS5', 'SBS6', 'SBS13', 'SBS15', 'SBS19', 'SBS30', 'SBS31', 'SBS32', 'SBS39', 'SBS44', 'SBS84', 'SBS86', 'SBS87', 'SBS96']
keep_columns_immucan = ['Type', 'SBS1', 'SBS2', 'SBS5', 'SBS13', 'SBS39', 'SBS87']
keep_columns_synergy = ['Type','SBS1', 'SBS2', 'SBS3', 'SBS5', 'SBS6', 'SBS13', 'SBS15', 'SBS19', 'SBS30', 'SBS31', 'SBS32', 'SBS39', 'SBS44', 'SBS84', 'SBS86', 'SBS87', 'SBS96']

df_immucan = df.loc[:, keep_columns_immucan]
df_synergy = df.loc[:, keep_columns_synergy]
print(len(df_immucan), len(df_synergy))

df_immucan.to_csv('IMMUcan_CUSTOM_COSMIC_v3.4_SBS_GRCh38_exome.txt', sep='\t', index=False)
df_synergy.to_csv('SYNERGY_CUSTOM_COSMIC_v3.4_SBS_GRCh38_exome.txt', sep='\t', index=False)


