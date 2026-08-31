from SigProfilerAssignment import Analyzer as Analyze

### IMMUcan before filtering - SBS

Analyze.cosmic_fit(samples="../data/counts/IMMUcan_BC1.SBS96.all",
                   output="./IMMUcan_before_filtering",
                   input_type="matrix",
                   context_type="96",
                   cosmic_version=3.4,
                   exome=True,
                   genome_build="GRCh38",
                   verbose=False)


### SYNERGY before filtering - SBS

Analyze.cosmic_fit(samples="../data/counts/Synergy_BC1.SBS96.all",
                   output="./SYNERGY_before_filtering",
                   input_type="matrix",
                   context_type="96",
                   cosmic_version=3.4,
                   exome=True,
                   genome_build="GRCh38",
                   verbose=False)