import pandas as pd


def add_mutation_counts(intervals_df, mutations_df):

    not_matched_df = pd.DataFrame({"chrom": [], 'position': [], 'gene': [], 'type': []})
    
    intervals_df['silent_count'] = 0
    intervals_df['missense_count'] = 0 

    mutations_index = 0
    mutations_len = len(mutations_df)

    for index, row in intervals_df.iterrows():

        while mutations_index < mutations_len and mutations_df.iloc[mutations_index]['Start_Position']-1 < row['start']:
            mutations_row = mutations_df.iloc[mutations_index]
            not_matched_df.loc[len(not_matched_df)] = [mutations_row['Chromosome'], mutations_row['Start_Position']-1, mutations_row['Hugo_Symbol'], mutations_row['Variant_Classification']]
            mutations_index += 1

        while mutations_index < mutations_len and \
            mutations_df.iloc[mutations_index]['Start_Position']-1 in range(int(row['start']), int(row['end'])+1):

            if mutations_index % 10000 == 0:
                print(mutations_index)

            if mutations_df.iloc[mutations_index]['Variant_Classification'] == 'Silent':
                intervals_df.at[index, 'silent_count'] += 1
            elif mutations_df.iloc[mutations_index]['Variant_Classification'] == 'Missense_Mutation':
                intervals_df.at[index, 'missense_count'] += 1

            mutations_index += 1

    return intervals_df, not_matched_df


def prep_mutations(mutations_df):

    mutations_df = mutations_df.sort_values(by=["Chromosome", "Start_Position"], ascending=[True, True]).reset_index(drop=True)

    mutations_df = mutations_df.drop(['file_path', 'project_id', 'End_Position', 'Strand', 'Variant_Type', 'CCDS'], axis=1)

    return mutations_df


def main():
    mutations_df = pd.read_csv("input_data/mutations.csv", sep="\t")
    regions_df = pd.read_csv("ref_data/modified/CCDS_split.csv")

    mutations_df = prep_mutations(mutations_df)

    mutations_df.to_csv("results/mutations.csv")

    # Initialize the DataFrames that will be output
    counts_df = regions_df.iloc[:0].copy()
    counts_df['silent_count'] = 0
    counts_df['missense_count'] = 0 
    print(counts_df)

    not_matched_df = pd.DataFrame({"chrom": [], 'position': [], 'gene': [], 'type': []})

    # Create a list of chromosomes in the CCDS DataFrame
    chroms = regions_df['chrom'].unique()

    for chrom in chroms:
        # Create temporary Dataframes consisting of only a single chromosme of data
        tmp_regions_df = regions_df[regions_df['chrom'] == chrom]
        tmp_mutations_df = mutations_df[mutations_df['Chromosome'] == chrom]

        # Run the split_chrom function to split the rows by SAE and SCE regions
        df1, df2 = add_mutation_counts(tmp_regions_df, tmp_mutations_df)

        # Concatenate regions to the output DataFrames
        counts_df = pd.concat([counts_df, df1])
        not_matched_df = pd.concat([not_matched_df, df2])
    
    # Return the Dataframes with reset indices
    counts_df = counts_df.reset_index(drop=True)
    not_matched_df = not_matched_df.reset_index(drop=True)


    counts_df.to_csv("results/counts.csv")
    not_matched_df.to_csv("results/not_matched.csv")


if __name__ == "__main__":
    main()