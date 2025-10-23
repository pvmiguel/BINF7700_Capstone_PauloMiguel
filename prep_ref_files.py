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

    return df


def CCDS_prep(path):
    # Load CCDS input file into dataframe
    df = pd.read_csv("ref_data/raw/CCDS.current.txt", sep="\t")

    df = df[df['ccds_status'] == 'Public']

    df["cds_locations"] = df["cds_locations"].str.replace("[", "")
    df["cds_locations"] = df["cds_locations"].str.replace("]", "")

    # Define columns to be exploded
    columns = ["cds_locations"]

    # Explode Columns
    df = explode_columns(df, columns)

    df[["cds_start", "cds_end"]] = df["cds_locations"].str.split('-', n=1, expand=True)

    df.insert(0, "chrom", "chr" + df["#chromosome"])

    df.drop(['#chromosome', 'cds_from', 'cds_to', 'match_type', 'ccds_status', 'cds_locations'], axis=1, inplace=True)

    return df



def main():
    SAE_df = SAE_SCE_prep("ref_data/raw/SAE.csv")
    SCE_df = SAE_SCE_prep("ref_data/raw/SCE.csv")
    CCDS_df = CCDS_prep("ref_data/raw/CCDS.current.txt")

    SAE_df.to_csv("ref_data/modified/SAE.csv", index=False)
    SCE_df.to_csv("ref_data/modified/SCE.csv", index=False)
    CCDS_df.to_csv("ref_data/modified/CCDS.csv", index=False)


if __name__ == "__main__":
    main()