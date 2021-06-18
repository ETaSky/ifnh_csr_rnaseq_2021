# Analysis of RNAseq data of IFNH CSR projects

## 0. Preparations
### Project folder structure
project
	|- README           # the top level description of content (this doc) and steps for analyzing the data
	|- CONTRIBUTING     # instructions for how to contribute to this project
	|- LICENSE          # the license for this project
    |
    |- code/            # scripts, notebooks used for analyze the data
    |
    |- data/            # raw and primary data, are not changed once created
    | |- raw/           # raw data and md5sum list, will not be altered
    | | |- fastqc/      # fastqc result for raw data
    | |- qc/            # quality filtered data by Trimmomatic
    | |- reference/     # reference genomes/annotations
    | |- mapped_STAR/   # mapping of qc data by STAR
    | |- mapped_salmon/ # mapping of qc data by salmon
    | |- featureCounts  # featureCounts output
    |
    |- doc/             # Documents that could be useful (e.g. Manuals)
    |
    |- logs/            # log files for analysis that would generate logs
    | |- 1/             # Parallel running trimmomatic sterr, stout for each files, only sterr is useful
