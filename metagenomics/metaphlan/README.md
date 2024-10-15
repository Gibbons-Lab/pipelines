## Workflow 

This pipeline takes single or paired-end fastq files as input and outputs taxonomic relative abundances


1. Adapter and quality trimming with fastp; quality reports in HTML and JSON are provided for each file
2. Taxonomy and abundance inference with Metaphlan4

---
 
## Environment set up and pipeline use

Install Anaconda

1. Download latest version of Anaconda: 
         
         wget https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh
         
2. Install using bash... for current AWS configuration Anaconda install dictory should be specified as /data/Anaconda3

         bash Anaconda3-2021.11-Linux-x86_64.sh

---
Pipeline is designed to run in a Anaconda environment and can be set up using the yaml file provided

1. Create ananconda environment with: 

        conda env create -f metaphlan_env.yml 
       
2. Set up directories for data and reference databases 
    Current set up assumes a format withe the following structure:   
    * A directory for the pipeline scripts (/data/pipeline)  
    * A directory for data and analysis (/data/analysis)  
     
     
3. Activate anaconda environment before use with:  

        conda activate metaphlan 

4. Download reference database (assumes you've activated the environment)
    		
		metaphlan --install         
        
5. Run pipeline with (see below for options):   

        nextflow run metaphlan.nf
  
---

## Options

    ~~~ Reference Based Metagenomics Workflow with Metaphlan4 ~~~

    Usage:
    A run using all,default parameters can be started with:
    > nextflow run ref_analyis.nf -resume

    A run with all parametrs set would look like:
    > nextflow run metaphlan_analysis.nf --data_dir=./data/analysis --fastq_dir=./data/fastq --single_end=false \\
                           --trim_front=5 --min_length=50 --quality_threshold=20 --read_length=150

    General options:
      --data_dir [str]              The main data directory for the analysis.
      --fastq_dir [str]             Directory containing raw data.
      --read_length [str]           The length of the reads.
      --single_end [bool]           Specifies that the input is single-end reads.
      --threads [int]               The maximum number of threads a single process can use.
                                    This is not the same as the maximum number of total threads used.
    Quality filter:
      --trim_front [str]            How many bases to trim from the 5' end of each read.
      --min_length [str]            Minimum accepted length for a read.
      --quality_threshold [str]     Smallest acceptable average quality.
