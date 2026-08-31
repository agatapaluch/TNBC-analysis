from SigProfilerAssignment import Analyzer as Analyze

### IMMUcan before filtering - DBS & ID

Analyze.cosmic_fit(samples="../data/counts/IMMUcan_BC1.DBS78.all",
                   output="DBS/DBS_IMMUcan_before_filtering",
                   input_type="matrix",
                   context_type="DBS78",
                   cosmic_version=3.4,
                   exome=True,
                   genome_build="GRCh38",
                   collapse_to_SBS96=False,
                   verbose=False)

Analyze.cosmic_fit(samples="../data/counts/IMMUcan_BC1.ID83.all",
                   output="ID/ID_IMMUcan_before_filtering",
                   input_type="matrix",
                   context_type="ID83",
                   cosmic_version=3.4,
                   exome=True,
                   genome_build="GRCh37", # ID data is available only for GRCh37
                   collapse_to_SBS96=False,
                   verbose=False)


### SYNERGY before filtering - DBS & ID

Analyze.cosmic_fit(samples="../data/counts/SYNERGY_BC1.DBS78.all",
                   output="DBS/DBS_SYNERGY_before_filtering",
                   input_type="matrix",
                   context_type="DBS78",
                   cosmic_version=3.4,
                   exome=True,
                   genome_build="GRCh38",
                   collapse_to_SBS96=False,
                   verbose=False)

Analyze.cosmic_fit(samples="../data/counts/SYNERGY_BC1.ID83.all",
                   output="ID/ID_SYNERGY_before_filtering",
                   input_type="matrix",
                   context_type="ID83",
                   cosmic_version=3.4,
                   exome=True,
                   genome_build="GRCh37", # ID data is available only for GRCh37
                   collapse_to_SBS96=False,
                   verbose=False)