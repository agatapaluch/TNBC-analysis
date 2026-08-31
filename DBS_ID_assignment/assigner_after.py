from SigProfilerAssignment import Analyzer as Analyze

### IMMUcan after filtering - DBS & ID

Analyze.cosmic_fit(samples="../data/counts/IMMUcan_BC1.DBS78.all",
                   output="DBS/DBS_IMMUcan_after_filtering",
                   input_type="matrix",
                   context_type="DBS78",
                   cosmic_version=3.4,
                   exome=True,
                   genome_build="GRCh38",
                   collapse_to_SBS96=False,
                   signature_database= 'IMMUcan_CUSTOM_COSMIC_v3.4_DBS_GRCh38_exome.txt',
                   verbose=False)

Analyze.cosmic_fit(samples="../data/counts/IMMUcan_BC1.ID83.all",
                   output="ID/ID_IMMUcan_after_filtering",
                   input_type="matrix",
                   context_type="ID83",
                   cosmic_version=3.4,
                   exome=True,
                   genome_build="GRCh37", # ID data is available only for GRCh37
                   signature_database='IMMUcan_CUSTOM_COSMIC_v3.4_ID_GRCh37.txt', # exome database not available
                   collapse_to_SBS96=False,
                   verbose=False)

### SYNERGY after filtering - DBS & ID

Analyze.cosmic_fit(samples="../data/counts/SYNERGY_BC1.DBS78.all",
                   output="DBS/DBS_SYNERGY_after_filtering",
                   input_type="matrix",
                   context_type="DBS78",
                   cosmic_version=3.4,
                   exome=True,
                   genome_build="GRCh38",
                   collapse_to_SBS96=False,
                   signature_database='SYNERGY_CUSTOM_COSMIC_v3.4_DBS_GRCh38_exome.txt',
                   verbose=False)

Analyze.cosmic_fit(samples="../data/counts/SYNERGY_BC1.ID83.all",
                   output="ID/ID_SYNERGY_after_filtering",
                   input_type="matrix",
                   context_type="ID83",
                   cosmic_version=3.4,
                   exome=True,
                   genome_build="GRCh37", # ID data is available only for GRCh37
                   signature_database= 'SYNERGY_CUSTOM_COSMIC_v3.4_ID_GRCh37.txt', # exome database not available
                   collapse_to_SBS96=False,
                   verbose=False)