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
    | |- qc/            # quality filtered data by Trimmomatic
    | |- reference/     # reference genomes/annotations
    | |- STAR_mapped/   # mapping of qc data by STAR
    |
    |- doc/             # Documents that could be useful (e.g. Manuals)
    |
    |- logs/            # log files for analysis that would generate logs
    | |- 1/             # Parallel running trimmomatic sterr, stout for each files, only sterr is useful

### Create the analysis conda environment and install packages
This section is constantly updating until the final analysis is completed along the process since new packages may be needed.
```sh
conda create --name RNAseqIFNH
conda activate RNAseqIFNH
conda install -c bioconda fastqc=0.11.9 trimmomatic=0.39 star=2.7.9a Subread=2.0.1
conda install -c conda-forge pigz=2.6
```

**Key Resources installed**
This section may be later moved to a requirement file
|Resource|Version|Note|
|:------|:------|:------|
|FastQC|0.11.9|QC|
|trimmomatic|0.39|QC|
|STAR|2.7.9.a|RNA aligner|
|Subread|2.0.1|the package contains featureCounts|
|pigz|2.6|Utility|
|Mouse Genome|GENCODE M27 (05.05.21)|Ref. genome, http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/GRCm39.primary_assembly.genome.fa.gz|
|Mouse Gene Annotation|GENCODE M27 (05.05.21)|GTF file, http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/gencode.vM27.primary_assembly.annotation.gtf.gz|

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

# Compress QCed files since too big (pigz needs to be installed)
cd data/qc/
pigz -p 24 *.fastq &
pigz -p 6 *.log &
cd ../../
```

### Align with STAR
#### Creating reference genome index (mouse genome)
The reference genome for the host animal (Swiss Webster strain in this project) is downloaded from GENCODE (https://www.gencodegenes.org/mouse/). Accoding to Section 2.2.1 of STAR manual, the GENCODE files marked with PRI (primary is downloaded)
```sh
# Genome
wget -P data/reference/ http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/GRCm39.primary_assembly.genome.fa.gz
pigz -d data/reference/gencode.vM27.primary_assembly.annotation.gtf.gz 

# Gene annotation (has to be matched with Genome)
wget -P data/reference/ http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/gencode.vM27.primary_assembly.annotation.gtf.gz
pigz -d data/reference/GRCm39.primary_assembly.genome.fa.gz 
(RNAseqIFNH) [jingchengw@annotate ifnh_csr_rnaseq_2021]$ 

# Generating genome indexes
mkdir -p data/reference/STAR_genome_ind
STAR --runThreadN 12 --runMode genomeGenerate --genomeDir data/reference/STAR_genome_ind --genomeFastaFiles data/reference/GRCm39.primary_assembly.genome.fa --sjdbGTFfile data/reference/gencode.vM27.primary_assembly.annotation.gtf --sjdbOverhang 149 &>logs/STAR_1_genome_index.log &
```

#### Mapping
```sh
# create file manifest
awk 'BEGIN{OFS="\t"}{print "data/qc/"$1"_O1_paired.fastq.gz","data/qc/"$1"_O2_paired.fastq.gz","ID:"$1}' data/raw_filenames.txt > data/STAR_filemanifest.tsv

# for STAR almost 20 minutes per sample
STAR --runThreadN 24 --genomeDir data/reference/STAR_genome_ind \
--readFilesManifest data/STAR_filemanifest.tsv \
--readFilesCommand pigz -dc \
--outSAMstrandField intronMotif \
--outFilterIntronMotifs RemoveNoncanonical \
--chimSegmentMin 30 \
--outSAMtype BAM Unsorted \
--outSAMunmapped Within \
--outSAMattributes NH HI AS nM NM MD jM jI MC ch RG \
--outFileNamePrefix data/mapped_STAR/ \
--limitBAMsortRAM 50000000000 &>logs/STAR_2_mapping.log & 

```

#### Generate feature table
Count gene using *featureCounts*, a subpackage of *Subread*. *featureCounts* may not work as expected if the bam file is not sorted by name for PE reads. # This step is pretty fast.
```sh
featureCounts -a data/reference/gencode.vM27.primary_assembly.annotation.gtf \
-o data/featureCounts/featureCounts_out \
-s 0 \
-p -C \
-T 24 \
--byReadGroup \
data/mapped_STAR/Aligned.out.bam &>logs/featureCounts.log &
```