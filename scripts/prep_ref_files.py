"""
File   : prep_ref_files.py
Created: 15-Nov-2025 by Paulo Miguel

Prepare the reference file that will be used for mutation counting of cancer 
MAF files. Creates a single CSV file with a single line per genomic interval
representing an SAE, SCE or normal exonic region of a human gene.
"""

import pandas as pd


def explode_columns(df, columns, islist=False):
    """
    Takes a list of columns in a Pandas DataFrame that includes a list
    of comma separated values and explodes the DataFrame creating a row
    per value in the list.

    Parameters:
    -----------
    df: Pandas DataFrame to be exploded
    columns: List of column names that will be exploded into rows
    islist: Optional parameter that indicates whether the columns to be
        exploded are a list.
    
    Returns:
    --------
    df: Pandas DataFrame with exploded columns
    """

    # Turn comma separated string into list for each column
    if islist == False:
        for column in columns:
            df = df.assign(**{column:df[column].str.split(',')})

    # Use explode to create a separate row per list entry
    df = df.explode(columns)

    return df


def SAE_SCE_prep(path):
    """
    Prepares SAE and SCE files for combination with the CCDS dataset.

    Parameters:
    -----------
    path: Path to SAE or SCE dataset that was downloaded from the UCSC 
        genome browser
    
    Returns:
    --------
    df: Pandas DataFrame with relevant SAE or SCE data
    """

    # Load SAE or SCE csv file into dataframe
    df = pd.read_csv(path)

    # Define columns to be exploded
    columns = ["blockSizes", "chromStarts"]

    # Explode Columns
    df = explode_columns(df, columns)

    # Cast blocksize and chromstart columns to integer type
    df['blockSize'] = df['blockSizes'].astype(int)
    df['chromStarts'] = df['chromStarts'].astype(int)

    # Drop unneeded columns
    df.drop(['score', 'thickStart', 'thickEnd', 'reserved', 'blockCount', 'blockSizes'], axis=1, inplace=True)

    # Calculate chromstart and chromend for exploded colunms
    df['chromStart'] = df['chromStart'] + df['chromStarts']
    # Adjust chromend value to be 0-based, see https://genome.ucsc.edu/FAQ/FAQtracks
    # for details (start is 0-based, end is 1-based by default)
    df['chromEnd'] = df['chromStart'] + df['blockSize'] - 1

    # Drop unneeded columns
    df.drop(['chromStarts'], axis=1, inplace=True)

    df = df.sort_values(by=['chrom', 'chromStart'])

    return df


def CCDS_prep(path):
    """
    Prepares CCDS file for combination with the SAE and SCE datasets.

    Parameters:
    -----------
    path: Path to CCDS dataset that was downloaded NCBI NIH 
    
    Returns:
    --------
    df: Pandas DataFrame with relevant SAE or SCE data
    """

    # Load CCDS input file into dataframe
    df = pd.read_csv("ref_data/raw/CCDS.current.txt", sep="\t")

    # Filter out CCDS data that has been withdrawn or is not complete 
    df = df[df['ccds_status'] == 'Public']
    df = df[df['match_type'] != 'Partial']

    # Prepare locations column for explosion
    df["cds_locations"] = df["cds_locations"].str.replace("[", "")
    df["cds_locations"] = df["cds_locations"].str.replace("]", "")

    # Define columns to be exploded
    columns = ["cds_locations"]

    # Explode Columns
    df = explode_columns(df, columns)

    # Create start and end columns based on the locations column
    df[["cds_start", "cds_end"]] = df["cds_locations"].str.split('-', n=1, expand=True)

    # Create chromosome column in chr# format
    df.insert(0, "chrom", "chr" + df["#chromosome"])

    # Drop unneeded columns
    df.drop(['#chromosome', 'cds_from', 'cds_to', 'match_type', 'ccds_status', 'cds_locations'], axis=1, inplace=True)

    # Cast start and end columns to integer type
    df['cds_start'] = df['cds_start'].str.strip().astype('int')
    df['cds_end'] = df['cds_end'].str.strip().astype('int')

    # Merge rows with exact same chromosome, start and end 
    df = df.groupby(['chrom','nc_accession','gene','gene_id','cds_strand','cds_start','cds_end'], as_index=False).agg({'ccds_id':','.join}).reset_index(drop=True)
    df['ccds_id'] = df['ccds_id'].apply(lambda x: f"[{x}]")

    # Sort by chromosome, start and end
    df = df.sort_values(by=['chrom', 'cds_start', 'cds_end'], ascending=[True, True, True]).reset_index(drop=True)

    # Merge rows with overalapping intervals 
    df = merge_overlapping(df, 'chrom', 'cds_start', 'cds_end')

    return df


def merge_overlapping(df, chrom_label, start_label, end_label):
    """
    Merges rows in a Pandas DataFrame that have overlapping intervals defined
    by their chromosome, start and end.

    Parameters:
    -----------
    df: Pandas DataFrame that contains interval overlaps
    chrom_label: String containing the name of the Pandas DataFrame column 
        where chromosome number is kept
    start_label: String containing the name of the Pandas DataFrame column 
        where the start of genomic interval is kept
    end_label: String containing the name of the Pandas DataFrame column 
        where the end of genomic interval is kept
    
    Returns:
    --------
    df: Pandas DataFrame with merged overlapping intervals 
    """

    # Define variables
    start = 0
    end = 0
    chrom = "chr0"
    keep_list = []
 
    for index, row in df.iterrows():
        #Check for when the chromosome changes to the next one
        if row[chrom_label] > chrom:
            # append to the keep list and reset variables
            keep_list.append([index-1, chrom, start, end])
            start = 0
            end = 0
            chrom = row[chrom_label]

        # Check if the start label is smaller than the current end    
        if row[start_label] <= end:
            # Check if the end label is larger than the current end
            if row[end_label] > end:
                # If True, update the end variable
                end = row[end_label]
        # Otherwise, append and reset the start and end variables
        else:
            keep_list.append([index-1, row[chrom_label], start, end])
            start = row[start_label]
            end = row[end_label]

    # Remove the first element of the keep list, since it contains all zeros
    keep_list.pop(0)

    # Append the final entry to the keep list
    keep_list.append([index, chrom, start, end])

    # Remove the first element from every chromosome, similar to previously
    real_keep_list = []
    for entry in keep_list:
        if entry[2] != 0:
            real_keep_list.append(entry)

    # Use zip and list comprehension to create individual lists for the index, 
    # chromosome, start and end
    unzipped_columns = zip(*real_keep_list)
    index, chrom, start, end = [list(col) for col in unzipped_columns]

    # Keep only the indices that were found previously
    df = df.iloc[index]

    # Reset start and end columns in the df DataFrame
    df[start_label] = start
    df[end_label] = end

    return df


def split_CCDS_by_coverage(CCDS_df, SAE_SCE_df):
    """
    Function for feeding chromosome restricted CCDS Dataframes and SAE_SCE DataFrames to 
    a function that splits the CCDS rows if they are SAE or SCE regions.
    
    Parameters:
    -----------
    CCDS_df: Pandas DataFrame with sorted CCDS information
    SAE_SCE_df : Pandas DataFrame with sorted SAE and SCE information
    
    Returns:
    --------
    keep_df = Pandas DataFrame with split CCDS rows by SAE and SCE
    del_df = Pandas DataFrame with SAE and SCE regions that did not fit in any CCDS rows
    """

    # Initialize the DataFrames that will be output
    keep_df = pd.DataFrame({"chrom": [], 'start': [], 'end': [], 'gene': [], 'type': []})
    del_df = pd.DataFrame({"chrom": [], 'start': [], 'end': [], 'gene': [], 'type': []})

    # Create a list of chromosomes in the CCDS DataFrame
    chroms = CCDS_df['chrom'].unique()

    for chrom in chroms:
        # Create temporary Dataframes consisting of only a single chromosme of data
        tmp_CCDS_df = CCDS_df[CCDS_df['chrom'] == chrom]
        tmp_SAE_SCE_df = SAE_SCE_df[SAE_SCE_df['chrom'] == chrom]

        # Run the split_chrom function to split the rows by SAE and SCE regions
        df1, df2 = _split_chrom(tmp_CCDS_df, tmp_SAE_SCE_df)

        # Concatenate regions to the output DataFrames
        keep_df = pd.concat([keep_df, df1])
        del_df = pd.concat([del_df, df2])
    
    # Return the Dataframes with reset indices
    keep_df = keep_df.reset_index(drop=True)
    del_df = del_df.reset_index(drop=True)
    
    return keep_df, del_df


def _split_chrom(CCDS_df, SAE_SCE_df):
    """
    Splits the rows in the CCDS DataFrame if they are SAE or SCE regions.
    
    Parameters:
    -----------
    CCDS_df: Pandas DataFrame with sorted CCDS information
    SAE_SCE_df : Pandas DataFrame with sorted SAE and SCE information
    
    Returns:
    --------
    df = Pandas DataFrame with split CCDS rows by SAE and SCE
    """

    # Duplicate the last row of the DataFrame
    last_row = SAE_SCE_df.iloc[[-1]]
    SAE_SCE_df = pd.concat([SAE_SCE_df, last_row])

    # Define metadata for output DataFrames
    keep_df = pd.DataFrame({"chrom": [], 'start': [], 'end': [], 'gene': [], 'type': []})
    del_df = pd.DataFrame({"chrom": [], 'start': [], 'end': [], 'gene': [], 'type': []})

    # Initialize index for iteration
    SAE_SCE_index = 0

    # Calculate the length of the SAE/SCE DataFrame
    SAE_SCE_len = len(SAE_SCE_df)

    for index, row in CCDS_df.iterrows():

        CCDS_start = row['cds_start']

        # Iterate if the SAE or SCE region starts before the CCDS region end and ensure that the 
        # index is not out of bounds
        while SAE_SCE_index < SAE_SCE_len and SAE_SCE_df.iloc[SAE_SCE_index]['chromStart'] <= row['cds_end']:

            # If the CCDS region starts after the current SAE SCE region, save the SAE SCE information
            # to the deleted DataFrame and reset the CCDS_start variable
            if CCDS_start > SAE_SCE_df.iloc[SAE_SCE_index]['chromStart']:
                del_df.loc[len(del_df)] = [row['chrom'], SAE_SCE_df.iloc[SAE_SCE_index]['chromStart'], CCDS_start-1, row['gene'], SAE_SCE_df.iloc[SAE_SCE_index]['type']]
                CCDS_start = SAE_SCE_df.iloc[SAE_SCE_index]['chromStart']

            # Otherwise, if the CCDS region start is not equal to the SAE SCE region start, add a line to the keep
            # DataFrame with normal type     
            elif CCDS_start != SAE_SCE_df.iloc[SAE_SCE_index]['chromStart']:
                keep_df.loc[len(keep_df)] = [row['chrom'], CCDS_start, SAE_SCE_df.iloc[SAE_SCE_index]['chromStart']-1, row['gene'], 'normal']
            
            # If the SAE SCE region ends before th CCDS region end, add a line to the keep DataFrame
            # and update the CCDS_start to the position after the end of the region that was just added
            if SAE_SCE_df.iloc[SAE_SCE_index]['chromEnd'] <= row['cds_end']:
                keep_df.loc[len(keep_df)] = [row['chrom'], SAE_SCE_df.iloc[SAE_SCE_index]['chromStart'], SAE_SCE_df.iloc[SAE_SCE_index]['chromEnd'], row['gene'], SAE_SCE_df.iloc[SAE_SCE_index]['type']]
                CCDS_start = SAE_SCE_df.iloc[SAE_SCE_index]['chromEnd'] + 1
            
            # Otherwise, add a line to the keep DataFrame for the region until the CCDS end and add 
            # a line to the delete DataFrame with the rest of the interval
            else: 
                keep_df.loc[len(keep_df)] = [row['chrom'], SAE_SCE_df.iloc[SAE_SCE_index]['chromStart'], row['cds_end'], row['gene'], SAE_SCE_df.iloc[SAE_SCE_index]['type']]
                del_df.loc[len(del_df)] = [row['chrom'], row['cds_end']+1, SAE_SCE_df.iloc[SAE_SCE_index]['chromEnd'], row['gene'], SAE_SCE_df.iloc[SAE_SCE_index]['type']]

            # If the index is 1 from the end, add a normal region making up the rest of the CCDS interval
            if SAE_SCE_index == SAE_SCE_len - 2:
                keep_df.loc[len(keep_df)] = [row['chrom'], SAE_SCE_df.iloc[SAE_SCE_index]['chromEnd']+1, row['cds_end'], row['gene'], 'normal']

            # Iterate the SAE SCE index and pull the relevant DataFrame row
            SAE_SCE_index += 1
            

        # If the algorithm did not enter the previous while loop, add the entire CCDS intveral as a 
        # normal region to the keep DataFrame
        if CCDS_start <= row['cds_end']:
            keep_df.loc[len(keep_df)] = [row['chrom'], CCDS_start, row['cds_end'], row['gene'], 'normal']
            
    return keep_df, del_df


def SAE_SCE_merge(SAE_df, SCE_df):
    """
    Merge and deduplicate SAE and SCE DataFrames
    
    Parameters:
    -----------
    SAE_df: Pandas DataFrame with SAE information
    SCE_df : Pandas DataFrame with SCE information
    
    Returns:
    --------
    df = Pandas DataFrame with merge and deduplicated SAE and SCE information
    """

    SAE_df['type'] = 'SAE'
    SCE_df['type'] = 'SCE'

    SAE_df = SAE_df.sort_values(by=['chrom', 'chromStart', 'chromEnd'], ascending=[True, True, True]).reset_index(drop=True)
    SCE_df = SCE_df.sort_values(by=['chrom', 'chromStart', 'chromEnd'], ascending=[True, True, True]).reset_index(drop=True)

    SAE_df = merge_overlapping(SAE_df, 'chrom', 'chromStart', 'chromEnd')
    SCE_df = merge_overlapping(SCE_df, 'chrom', 'chromStart', 'chromEnd')

    df = pd.concat([SAE_df, SCE_df])

    df = df.sort_values(by=['chrom', 'chromStart', 'chromEnd'], ascending=[True, True, True]).reset_index(drop=True)

    return df


def main():
    SAE_df = SAE_SCE_prep("ref_data/raw/SAE.csv")
    SCE_df = SAE_SCE_prep("ref_data/raw/SCE.csv")
    CCDS_df = CCDS_prep("ref_data/raw/CCDS.current.txt")

    SAE_SCE_df = SAE_SCE_merge(SAE_df, SCE_df)

    CCDS_split_df, CCDS_del_df = split_CCDS_by_coverage(CCDS_df, SAE_SCE_df)

    SAE_df.to_csv("ref_data/modified/SAE.csv", index=False)
    SCE_df.to_csv("ref_data/modified/SCE.csv", index=False)
    SAE_SCE_df.to_csv("ref_data/modified/SAE_SCE.csv", index=False)
    CCDS_df.to_csv("ref_data/modified/CCDS.csv", index=False)
    CCDS_split_df.to_csv("ref_data/modified/CCDS_split.csv", index=False)
    CCDS_del_df.to_csv("ref_data/modified/CCDS_del.csv", index=False)


if __name__ == "__main__":
    main()