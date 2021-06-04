# Analysis of RNAseq data of IFNH CSR projects

## 0. Preparations
### Project folder structure
project
	|- README           # the top level description of content (this doc) and steps for analyzing the data
	|- CONTRIBUTING     # instructions for how to contribute to this project
	|- LICENSE          # the license for this project
    |
    |- data/            # raw and primary data, are not changed once created
    | |- raw/           # raw data and md5sum list, will not be altered
    | | |- fastqc/      # fastqc result for raw data
    | |- qc/            # Trimmomatic outputs
    |
    |- logs/            # log files for analysis that would generate logs
    | |- 1/             # Parallel running trimmomatic sterr, stout for each files, only sterr is useful

### Create the analysis conda environment and install packages
This section is constantly updating until the final analysis is completed along the process since new packages may be needed.
```sh
conda create --name RNAseqIFNH
conda activate RNAseqIFNH
conda install -c bioconda fastqc=0.11.9 trimmomatic=0.39 star=2.7.9a
```

**Packages/Softwares installed**
This section may be later moved to a requirement file
|Package|Version|
|:------|:------|
|FastQC|0.11.9|
|trimmomatic|0.39|
|STAR|2.7.9.a|

## 1. Quality checking and processing

### File integrity
Move raw data into project folder under `data/raw/`, generate a md5 list for all the files in the folder, and check with the md5sum_list.txt.

```sh
cd data/raw/
md5sum *.fastq.gz
```
All files currently in the `raw/` folder has been checked out with the md5sum list.

### Quality checking with fastqc
Run `fastqc`
```sh
fastqc -o data/raw/fastqc *.fastq.gz --noextract -t 16 &>logs/fastqc.logs &
```

### Quality trimming with Trimmomatic
```sh
# Generate a file name list
for f in data/raw/*R1_001.fastq.gz; do echo $f; done | awk -F "/" '{print $3}' | sed -e 's/_R1_001.fastq.gz//g' > data/raw_filenames.txt
# Using GNU parallel to perform trimmomatic
parallel -a data/raw_filenames.txt -j 3 --results logs/ --joblog logs/trimmomatic.log \
    trimmomatic PE -threads 8 \
        -phred33 \
        -trimlog data/qc/{}.log \
        data/raw/{}_R1_001.fastq.gz data/raw/{}_R2_001.fastq.gz \
        data/qc/{}_O1_paired.fastq data/qc/{}_O1_unpaired.fastq \
        data/qc/{}_O2_paired.fastq data/qc/{}_O2_unpaired.fastq \
        ILLUMINACLIP:$HOME/miniconda3/envs/RNAseqIFNH/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads \
        LEADING:20 TRAILING:20 \
        SLIDINGWINDOW:4:15 \
        MINLEN:35 &

# Check log files for length distribution
while IFS='' read -r line || [[ -n "$line" ]]; do \
    awk -v a="$line"  -F" " '(NR%2==1){F[$3]+=1}(NR%2==0){R[$3]+=1}END{for (i in F) print a,"R1",i,F[i]}END{for (j in R) print a,"R2",j,R[j]}' data/qc/${line}.log; \
done < data/raw_filenames.txt &>data/qc/length_distribution.txt &

# Check log files for summary quality trimming results
while IFS='' read -r line || [[ -n "$line" ]]; do \
    echo $line; grep "Input Read Pairs" logs/1/${line}/stderr;
done < data/raw_filenames.txt &>data/qc/summary_stats.txt &
```