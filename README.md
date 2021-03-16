# RNA2CM
RNA2CM is a tool for the identification of cancer-related mutations from RNA-seq data

## Installation
1| Nextflow and Docker (Singularity for execution on SLURM-clusters) are the only pre-requisites for the RNA2CM tool. Install both if needed and make sure they are properly running on your system. If the following commands do not generate any error message you are good to go.
```bash
nextflow run hello # test that nextflow is working
docker run hello-world # test that docker is working
```
2| Download and extract the project directory (using either git or wget):
```bash
git clone https://github.com/elyadlezmi/RNA2CM.git # clone the project using git
```
or
```bash
wget https://github.com/elyadlezmi/RNA2CM/archive/master.zip && unzip master.zip && mv RNA2CM-master RNA2CM  
```
3| Download the files CosmicMutantExportCensus.tsv.gz and CosmicCodingMuts.vcf.gz from the COSMIC website (login required), then move them into the project’s subdirectory named data (RNA2CM/data).

4| Execute the script named RNA2CMsetup.nf which is responsible for setting up all the reference data. This script can take two optional parameters:

--cpu: The number of threads for multi-threading (int, default 8).

--readLength: The expected Illumina read length for optimal alignment by STAR (int, default 100).

-profile: Choose the executor profile between a standard dockerized usage on a local workstation, usage on a SLURM cluster (requires Singularity instead of Docker) or a fully local execution which is the least recommended option (standard/cluster/local, default: standard).

Example for running the setup using 4 CPU and with a read length of 75bp:
```bash
RNA2CMsetup.nf --cpu 4 --readLength 75
```
This and the following examples assume that both the RNA2CM directory, and the Nextflow executable are within your PATH variable. If not, you should specify the nextflow interperter and state the absolute paths in your commands. e.g. if nextflow is located in /Home/apps and the project is located in /Home/bioinformatics, type:
```bash
/Home/apps/nextflow /Home/bioinformatics/RNA2CM/RNA2CMsetup.nf --cpu 4 --readLength 75
```

## Usage
```bash
RNA2CM.nf --fastq your_sample.fastq.gz # for single-end reads
RNA2CM.nf --fastq your_sample_1.fastq.gz --fastq2 your_sample_2.fastq.gz # for paired-ends reads
```
Optional arguments (Note that the only required arguments are RNA-seq reads, output is generated into the working directory):

--cpu:(int, default 8).

--keepInter: Whether to keep intermediate alignment and VCF files (true/false, default: false) 

--filterMouse: Whether to perform mouse contamination cleanup (true/false, default true)

-profile: (standard/cluster/local, default: standard).

Example for a paired-ends RNA-seq run, using 4 CPUs, keeping intermediate files:
```bash
RNA2CM.nf --fastq esc_1.fastq.gz --fastq2 esc_2.fastq.gz --cpu 4 --keepInter true 
```
Example for a single-ends RNA-seq run, skipping mouse read filtration and running on a SLURM cluster:
```bash
RNA2CM.nf --fastq SRR3090631.fastq.gz --filterMouse false -profile cluster
```
