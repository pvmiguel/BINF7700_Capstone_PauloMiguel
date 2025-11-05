import argparse
import pandas as pd
from genome_kit import Genome
from genome_kit import Interval


def get_cli_args():
    """
    Get command line inputs for analysis

    Parameters:
    None

    Returns:
    parser.parse_args(): Instance of argparse arguments
    """

    # Create parser object
    parser = argparse.ArgumentParser()

    # Add infile argument to parser
    parser.add_argument('-i',
                        '--infile',
                        type=str,
                        help='Path to the cancer MAF file to open',
                        default="input_data/229a0a22-895b-4a65-bd13-3532c7874519.wxs.aliquot_ensemble_masked.maf")

    return parser.parse_args()


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
    
    df = df[df['Variant_Classification'].isin(['Missense_Mutation', 'Silent'])]

    return df


def check_nucleotide(genome, chrom, pos, strand, ref_allele):

    interval = Interval(chrom, strand, pos-1, pos, 'hg38')

    print(genome.dna(interval), ref_allele)

    return genome.dna(interval) == ref_allele

def main():

    # Get arguments from the command line and set variables
    args = get_cli_args()
    infile = args.infile

    # Load cancer MAF file
    sample_df = import_maf(infile)

    genome = Genome('hg38')

    for index, row in sample_df.iterrows():
        if not check_nucleotide(genome, row['Chromosome'], row['Start_Position'], row['Strand'], row['Reference_Allele']):
            print(row)


if __name__ == "__main__":
    main()