"""
File   : validate_mutations.py
Created: 15-Nov-2025 by Paulo Miguel

Validate mutations from MAF files and output them into a single CSV file.
"""

import pandas as pd
import logging
import sys
import csv
import os
from genome_kit import Genome
from genome_kit import Interval
from datetime import datetime


def import_maf(path):
    """
    Import MAF file
    
    Parameters:
    -----------
    path: Path to the MAF file to import
    
    Returns:
    --------
    df: Pandas Dataframe containing MAF file information
    """

    # Read the MAF file, ignoring the first 7 rows
    df = pd.read_csv(path, sep='\t', header=7)

    # Keep only relevant columns
    df = df[['Hugo_Symbol',
            'Entrez_Gene_Id',
            'NCBI_Build',
            'Chromosome',
            'Start_Position',
            'End_Position',
            'Strand',
            'Variant_Classification',
            'Variant_Type',
            'Reference_Allele',
            'Tumor_Seq_Allele1',
            'Tumor_Seq_Allele2',
            'Mutation_Status',
            'Gene',
            'Feature',
            'cDNA_position',
            'CDS_position',
            'Protein_position',
            'Amino_acids',
            'Codons',
            'CCDS']]
    
    return df


def import_metadata(path):
    """
    Import MAF metadata csv file
    
    Parameters:
    -----------
    path: Path to the MAF metadata CSV file to import
    
    Returns:
    --------
    metadata: list containing MAF metadata
    """

    # Try opening metadata
    try:
        with open(path, 'r') as file:

            # Create list of metadata data
            csv_reader = csv.reader(file)
            metadata = list(csv_reader)
            # Remove the header
            metadata.pop(0)

        return metadata
    
    # If FileNotFound, return error
    except FileNotFoundError:
        print(f"Error: {path} not found. Please ensure the file exists.")


def check_nucleotide(genome, chrom, pos, strand, ref_allele):
    """
    Check if a nucleotide matches the reference nucleotide for a give genomic position
    
    Parameters:
    -----------
    genome: genome kit Genome object for the genome in question
    chrom: Chromosome
    pos: Location of the nucleotide within the chromosome to be checked
    strand: The direction of the genome the nucleotide came from
    ref_allele: The nucleotide that is being checked
    
    Returns:
    --------
    Boolean: True if there is a match, False if it does not match
    """

    # Find the nucleotide at the genomic location
    interval = Interval(chrom, strand, pos-1, pos, 'hg38')

    # Return True or False based on if the nucletides match
    return genome.dna(interval) == ref_allele


def main():
    """
    Main Process Flow
    """

    # Set up logging
    logfile = "mutationvalidation.log"

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s  %(message)s',
        handlers=[
            logging.FileHandler(f"logs/{datetime.now().strftime("%Y%m%d")}_{logfile}"),
            logging.StreamHandler(sys.stdout)
        ]
    )

    # Create filehandle for the mutations CSV file and add the header
    fh_mutations = open("input_data/mutations.csv", "w")
    fh_mutations.write("file_path\tproject_id\tproject_name\tdisease_type\tprimary_site\tHugo_Symbol\tEntrez_Gene_Id\tChromosome\tStart_Position\tEnd_Position\tStrand\tVariant_Classification\tVariant_Type\tCCDS\n")

    # Create filehandle for the raw counts CSV file and add the header
    fh_counts = open("input_data/raw_counts.csv", "w")
    fh_counts.write("file_path\tsaved_mutations\ttotal_mutations\n")

    logger = logging.getLogger(__name__) 
    
    # Import the MAF metadata
    metadata = import_metadata("maf_files/maf_metadata.csv")

    # Initialize mutation counters
    total_mutations = 0
    saved_mutations = 0

    # Iterate through the files found in the metadata file
    for line in metadata:

        # Import the MAF file
        sample_df = import_maf(f"maf_files/{line[0]}")

        # Initialize the genome kit Genome object
        genome = Genome('hg38')

        # Initialize the counts for the file
        file_total_mutations = 0
        file_saved_mutations = 0

        logger.info(f"{'='*60}")
        logger.info(f"File: {line[0]}")

        # If the file is whole exon sequencing
        if line[4] == "WXS":

            # Iterate through the rows of the MAF
            for index, row in sample_df.iterrows():
                
                # Reject and log mutations that are not missense or silent SNPs
                if row['Variant_Classification'] not in ['Missense_Mutation', 'Silent'] or row['Variant_Type'] != 'SNP':
                    logger.info(f"    Index: {index} not a Silent or Missense SNP")

                # Check if the allele nucleotide matches the reference nucleotide
                elif check_nucleotide(genome, row['Chromosome'], row['Start_Position'], row['Strand'], row['Reference_Allele']):

                    # Write the mutation to the mutations file
                    fh_mutations.write(f"{line[0]}\t{line[6]}\t{line[7]}\t{line[8]}\t{line[9]}\t{row['Hugo_Symbol']}\t{row['Entrez_Gene_Id']}\t{row['Chromosome']}\t{row['Start_Position']}\t{row['End_Position']}\t{row['Strand']}\t{row['Variant_Classification']}\t{row['Variant_Type']}\t{row['CCDS']}\n")
                    
                    # Add one the file saved mutations
                    file_saved_mutations += 1

                # If the allele nucleotide does not match, do not add
                else:
                    logger.info(f"    Index: {index} DOES NOT MATCH: Reference Allele ({row['Reference_Allele']}) does not match Genome hg38")

                # Add one to the file total mutations
                file_total_mutations += 1

            # Log number of mutations saved and total 
            logger.info(f"File {line[0]} Saved Mutations = {file_saved_mutations}")
            logger.info(f"File {line[0]} Total Mutations = {file_total_mutations}")
            logger.info(f"{'='*60}\n")

            # Write to the raw count file the number of saved and total mutations for the file
            fh_counts.write(f"{line[0]}\t{file_saved_mutations}\t{file_total_mutations}\n")

            # Add the file total and saved to the total and saved variables and reset the file variables
            saved_mutations += file_saved_mutations
            total_mutations += file_total_mutations
            file_saved_mutations = 0
            file_total_mutations = 0
        
        # If not WXS, skip
        else: 
            logger.info(f"    File {line[0]} is not a Whole Exome Sequencing (WXS) project and was not saved")
            logger.info(f"{'='*60}\n")
    
    # Write number of saved and total mutations across all files
    logger.info(f"Saved {saved_mutations} mutations out of {total_mutations}")

    # Close filehandles
    fh_mutations.close()
    fh_counts.close()


if __name__ == "__main__":
    main()