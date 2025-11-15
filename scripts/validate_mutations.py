import pandas as pd
import logging
import sys
import csv
import os
from genome_kit import Genome
from genome_kit import Interval
from datetime import datetime


def import_maf(path):

    df = pd.read_csv(path, sep='\t', header=7)

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

    try:
        with open(path, 'r') as file:
            csv_reader = csv.reader(file)
            metadata = list(csv_reader)
            metadata.pop(0)
        return metadata
    except FileNotFoundError:
        print(f"Error: {path} not found. Please ensure the file exists.")


def check_nucleotide(genome, chrom, pos, strand, ref_allele):

    interval = Interval(chrom, strand, pos-1, pos, 'hg38')

    return genome.dna(interval) == ref_allele

def main():

    logfile = "mutationvalidation.log"

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s  %(message)s',
        handlers=[
            logging.FileHandler(f"logs/{datetime.now().strftime("%Y%m%d")}_{logfile}"),
            logging.StreamHandler(sys.stdout)
        ]
    )

    fh_mutations = open("input_data/mutations.csv", "w")
    fh_mutations.write("file_path\tproject_id\tproject_name\tdisease_type\tprimary_site\tHugo_Symbol\tEntrez_Gene_Id\tChromosome\tStart_Position\tEnd_Position\tStrand\tVariant_Classification\tVariant_Type\tCCDS\n")

    fh_counts = open("input_data/raw_counts.csv", "w")
    fh_counts.write("file_path\tsaved_mutations\ttotal_mutations\n")

    logger = logging.getLogger(__name__) 
    
    metadata = import_metadata("maf_files/maf_metadata.csv")

    total_mutations = 0
    saved_mutations = 0

    for line in metadata:

        sample_df = import_maf(f"maf_files/{line[0]}")

        genome = Genome('hg38')

        file_total_mutations = 0
        file_saved_mutations = 0

        logger.info(f"{'='*60}")
        logger.info(f"File: {line[0]}")

        if line[4] == "WXS":

            for index, row in sample_df.iterrows():
                
                if row['Variant_Classification'] not in ['Missense_Mutation', 'Silent'] or row['Variant_Type'] != 'SNP':
                    logger.info(f"    Index: {index} not a Silent or Missense SNP")

                elif check_nucleotide(genome, row['Chromosome'], row['Start_Position'], row['Strand'], row['Reference_Allele']):
                    fh_mutations.write(f"{line[0]}\t{line[6]}\t{line[7]}\t{line[8]}\t{line[9]}\t{row['Hugo_Symbol']}\t{row['Entrez_Gene_Id']}\t{row['Chromosome']}\t{row['Start_Position']}\t{row['End_Position']}\t{row['Strand']}\t{row['Variant_Classification']}\t{row['Variant_Type']}\t{row['CCDS']}\n")
                    
                    file_saved_mutations += 1

                else:
                    logger.info(f"    Index: {index} DOES NOT MATCH: Reference Allele ({row['Reference_Allele']}) does not match Genome hg38")

                file_total_mutations += 1

            logger.info(f"File {line[0]} Saved Mutations = {file_saved_mutations}")
            logger.info(f"File {line[0]} Total Mutations = {file_total_mutations}")
            logger.info(f"{'='*60}\n")

            fh_counts.write(f"{line[0]}\t{file_saved_mutations}\t{file_total_mutations}\n")

            saved_mutations += file_saved_mutations
            total_mutations += file_total_mutations
            file_saved_mutations = 0
            file_total_mutations = 0
        
        else: 
            logger.info(f"    File {line[0]} is not a Whole Exome Sequencing (WXS) project and was not saved")
            logger.info(f"{'='*60}\n")
    
    logger.info(f"Saved {saved_mutations} mutations out of {total_mutations}")

    fh_mutations.close()
    fh_counts.close()


if __name__ == "__main__":
    main()