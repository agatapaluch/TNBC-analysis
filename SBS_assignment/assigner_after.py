from SigProfilerAssignment import Analyzer as Analyze

### IMMUcan after filtering - SBS

Analyze.cosmic_fit(samples="../data/counts/IMMUcan_BC1.SBS96.all",
                   output="./IMMUcan_after_filtering",
                   input_type="matrix",
                   context_type="96",
                   cosmic_version=3.4,
                   exome=True,
                   signature_database = 'IMMUcan_CUSTOM_COSMIC_v3.4_SBS_GRCh38_exome.txt',
                   genome_build="GRCh38",
                   verbose=False)

### SYNERGY after filtering - SBS

Analyze.cosmic_fit(samples="../data/counts/Synergy_BC1.SBS96.all",
                   output="./SYNERGY_after_filtering",
                   input_type="matrix",
                   context_type="96",
                   cosmic_version=3.4,
                   exome=True,
                   genome_build="GRCh38",
                   signature_database = 'SYNERGY_CUSTOM_COSMIC_v3.4_SBS_GRCh38_exome.txt',
                   verbose=False)