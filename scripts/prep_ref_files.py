import pandas as pd


def explode_columns(df, columns, islist=False):
    # Turn comma seprated string into list for each column
    if islist == False:
        for column in columns:
            df = df.assign(**{column:df[column].str.split(',')})

    # Use explode to create a separate row per list entry
    df = df.explode(columns)

    return df


def SAE_SCE_prep(path):
    # Load SAE or SCE csv file into dataframe
    df = pd.read_csv(path)

    # Define columns to be exploded
    columns = ["blockSizes", "chromStarts"]

    # Explode Columns
    df = explode_columns(df, columns)

    df['blockSize'] = df['blockSizes'].astype(int)
    df['chromStarts'] = df['chromStarts'].astype(int)

    df.drop(['score', 'thickStart', 'thickEnd', 'reserved', 'blockCount', 'blockSizes'], axis=1, inplace=True)

    df['chromStart'] = df['chromStart'] + df['chromStarts']
    df['chromEnd'] = df['chromStart'] + df['blockSize']

    df.drop(['chromStarts'], axis=1, inplace=True)

    df = df.sort_values(by=['chrom', 'chromStart'])

    return df


def CCDS_prep(path):
    # Load CCDS input file into dataframe
    df = pd.read_csv("ref_data/raw/CCDS.current.txt", sep="\t")

    df = df[df['ccds_status'] == 'Public']
    df = df[df['match_type'] != 'Partial']

    df["cds_locations"] = df["cds_locations"].str.replace("[", "")
    df["cds_locations"] = df["cds_locations"].str.replace("]", "")

    # Define columns to be exploded
    columns = ["cds_locations"]

    # Explode Columns
    df = explode_columns(df, columns)

    df[["cds_start", "cds_end"]] = df["cds_locations"].str.split('-', n=1, expand=True)

    df.insert(0, "chrom", "chr" + df["#chromosome"])

    df.drop(['#chromosome', 'cds_from', 'cds_to', 'match_type', 'ccds_status', 'cds_locations'], axis=1, inplace=True)

    df['cds_start'] = df['cds_start'].str.strip().astype('int')
    df['cds_end'] = df['cds_end'].str.strip().astype('int')

    df = df.groupby(['chrom','nc_accession','gene','gene_id','cds_strand','cds_start','cds_end'], as_index=False).agg({'ccds_id':','.join}).reset_index(drop=True)

    df['ccds_id'] = df['ccds_id'].apply(lambda x: f"[{x}]")

    df = df.sort_values(by=['chrom', 'cds_start', 'cds_end'], ascending=[True, True, True]).reset_index(drop=True)

    df = merge_overlapping_CCDS(df, 'chrom', 'cds_start', 'cds_end')

    return df


def merge_overlapping_CCDS(df, chrom_label, start_label, end_label):
    start = 0
    end = 0
    chrom = "chr0"
    keep_list = []
    for index, row in df.iterrows():
        if row[chrom_label] > chrom:
            keep_list.append([index-1, chrom, start, end])
            start = 0
            end = 0
            chrom = row[chrom_label]
        if row[start_label] <= end:
            if row[end_label] > end:
                end = row[end_label]
        else:
            keep_list.append([index-1, row[chrom_label], start, end])
            start = row[start_label]
            end = row[end_label]

    keep_list.pop(0)
    keep_list.append([index, chrom, start, end])

    real_keep_list = []
    for entry in keep_list:
        if entry[2] != 0:
            real_keep_list.append(entry)

    unzipped_columns = zip(*real_keep_list)

    index, chrom, start, end = [list(col) for col in unzipped_columns]

    df = df.iloc[index]
    df[start_label] = start
    df[end_label] = end

    return df


def CCDS_split(CCDS_df, SAE_df, SCE_df):

    SAE_df['Type'] = 'SAE'
    SCE_df['Type'] = 'SCE'

    SAE_SCE_df = pd.concat([SAE_df, SCE_df])

    return True


def main():
    SAE_df = SAE_SCE_prep("ref_data/raw/SAE.csv")
    SCE_df = SAE_SCE_prep("ref_data/raw/SCE.csv")
    CCDS_df = CCDS_prep("ref_data/raw/CCDS.current.txt")

    SAE_df.to_csv("ref_data/modified/SAE.csv", index=False)
    SCE_df.to_csv("ref_data/modified/SCE.csv", index=False)
    CCDS_df.to_csv("ref_data/modified/CCDS.csv", index=False)


if __name__ == "__main__":
    main()