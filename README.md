# BINF7700 Capstone Project Repository
### Author: Paulo Miguel

## Purpose

This repository contains an analysis of cancer mutations with respect to synonymous contraint elements (SCEs) and sysnonymous acceleration elements (SAEs). It contains all the code that was used to transform raw CCDS, SAE and SCE data into genomic intervals and the counting of mutations from the GDC Data Portal into those regions. A statistical analysis and other data exploration is conducted in a series of Jupyter notebooks. 

## Dependencies

All environment dependencies for this analysis can be found in the environment.yml file in the repository. The dependencies are as follows:

* python=3.12.1
* numpy=2.3.5
* pandas=2.3.3
* genomekit=7.2.2
* sh=2.2.2
* jupyter=1.1.1
* matplotlib=3.10.8
* seaborn=0.13.2
* scipy=1.16.3
* statsmodels=0.14.5
* gseapy=1.1.11
* adjusttext=1.3.0
* dna_features_viewer

Installation of the dependencies using the environment.yml file requires conda. Run the below command to create the conda environment.

```
conda env create -f environment.yml
```

## Data

### Availability of Raw Reference Files

NOTE: The raw files used in this analysis have been provided in this directory. They can be found zipped in `ref_data/raw.zip`.

The location and date of download of all the reference files used in this analysis can be found below:

Consensus CDS (CCDS) file: 
* Downloaded on 2025/10/22 from https://ftp.ncbi.nih.gov/pub/CCDS/current_human/
* CCDS.current.txt file was downloaded
* Date of last modification of file on download date was 2025/09/28 at 13:12

Synonymous Constraint Element file:
* Downloaded on 2025/10/16 from https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=hub_377623&hgta_track=hub_377623_SCE&hgta_table=hub_377623_SCE&hgta_doSchema=describe+table+schema

Synonymous Acceleration Element file:
* Downloaded on 2025/10/16 from https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=hub_377623&hgta_track=hub_377623_SAE&hgta_table=hub_377623_SAE&hgta_doSchema=describe+table+schema

### Modification of Raw Reference Files

NOTE: The modified files used in this analysis have been provided in this directory. They can be found zipped in `ref_data/modified.zip`.

The modified files used in this analysis were generated using the `prep_ref_files.py` script which can be found in the `scripts/` folder of this repository. To run the scipt, go to root file of the repository and run `prep_ref_files.py`.

```
python scripts/prep_ref_files.py
```

### Downloading Mutation Data

The mutation data used in this analysis were sourced from the National Cancer Institute GDC Data Portal (https://portal.gdc.cancer.gov/). To download the files locally, the `maf_downloader.py` script was created. It can be ran using the following command.

```
python scripts/count_mutations.py -n [number of files] -l [logfile name]
```

NOTE: This script has not yet been modified to recognize which files have already been downloaded to the repository and download other undownloaded files. To get all open access WXS MAF files on the GDC Data Portal, use -n 20000. Ideally, use screen or sbatch, since the download will take 8+ hours.

### Mutation Validation

All mutation data downloaded used in this analysis was verified, by checking if the reference allele in the files matched HG38, and filtered, keeping only the silent (synonymous) and missense (non-synonymous) single nucleotide polymorphisms (SNP). Following the above step of downloading the mutation files, this step can be completed using the `validate_mutations.py` script. It can be ran using the following command.

```
python scripts/validate_mutations.py
```

### Mutation Counting

The validation step returns a CSV file with all validated mutations that can then be counted using the regions created using the `prep_ref_files.py` script. The mutation counting was done using the `count_mutations.py` script and can be ran using the following command.

```
python scripts/count_mutations.py
```

## Analysis

The notebooks folder includes several Jupyter notebooks used for various analyses. They are are follows:

### analysis.ipynb

This Jupyter notebook contains exploratory analyses into the count data, hypothesis testing on the mutation distributions in the genes and over-representation analysis (ORA) on the statistically significant genes.

### gene_distributions.ipynb

This Jupyter notebook contains visualizations of the mutation rate over the length of individual genes. It is used to visualize and determine qualitatively which regions of a gene have higher or lower mutation rates.

### not_matched.ipynb

This Jupyter notebook contains an analysis on the mutations that did not match the CCDS regions from the raw files. 

### sample_dist.ipynb

This Jupyter notebook contains the creating of a visualization for the number of mutations found in each of the analyzed samples.

## Acknowledgements

This project was completed as part of the capstone project course (BINF 7700) for the Master of Bioinformatics program at Northeastern University Toronto. 

Special thank you to:
* Max Wolf (Supervisor)
* Oyeronke Ayansola (Course Instructor)

## License

See the [LICENSE](LICENSE) file