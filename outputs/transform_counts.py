import gzip
import os
import pandas as pd
import re


def strip_version(ensembl_id):
    return ensembl_id.split('.')[0]


def get_gene_lengths_from_gtf(filepath):
    """
    Parse a GTF file and return a DataFrame with gene_id and their exon-based lengths.

    Parameters:
    filepath (str): Path to the GTF file (.gtf or .gtf.gz)

    Returns:
    pd.DataFrame: DataFrame with columns ['gene_id', 'gene_length']
    """
    def parse_attributes(attr_str):
        attrs = {}
        for attr in attr_str.strip().split(';'):
            if attr.strip() == '':
                continue
            key, value = attr.strip().split(' ', 1)
            attrs[key] = value.strip('"')
        return attrs

    exon_intervals = {}

    open_func = gzip.open if filepath.endswith('.gz') else open
    with open_func(filepath, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) != 9:
                continue
            feature = parts[2]
            if feature != 'exon':
                continue
            start, end = int(parts[3]), int(parts[4])
            attributes = parse_attributes(parts[8])
            gene_id = attributes.get('gene_id')
            if gene_id is None:
                continue
            if gene_id not in exon_intervals:
                exon_intervals[gene_id] = []
            exon_intervals[gene_id].append((start, end))

    # Merge overlapping exon intervals
    gene_lengths = []
    for gene_id, intervals in exon_intervals.items():
        intervals.sort()
        merged = []
        for start, end in intervals:
            if not merged or merged[-1][1] < start:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)
        total_length = sum(end - start + 1 for start, end in merged)
        gene_lengths.append((gene_id, total_length))

    df = pd.DataFrame(gene_lengths, columns=['gene_id', 'gene_length'])
    return df


def process_featureCounts(filepath):
    """Process featureCounts output file."""
    print(f"\t Processing featureCounts file: {filepath}")
    df = pd.read_csv(filepath, sep="\t", comment="#")
    df = df.iloc[:, [0, -2, -1]]  # Keep gene ID, length, and count columns
    df.columns = ["gene_id", "length", "count"]
    df["gene_id"] = df["gene_id"].apply(strip_version)  # Strip version from gene_id
    df["TPM"] = (df["count"] / df["length"]) / (df["count"] / df["length"]).sum() * 1e6  # Calculate TPM
    df = df[["gene_id", "TPM"]]  # Keep only gene_id and TPM
    return df

def process_htseq(filepath, gene_lengths_df):
    """
    Process HTSeq output and compute TPM values using gene lengths.

    Parameters:
    filepath (str): Path to HTSeq count file (TSV)
    gene_lengths_df (pd.DataFrame): DataFrame with 'gene_id' and 'gene_length' (in bp)

    Returns:
    pd.DataFrame: DataFrame with ['gene_id', 'counts', 'gene_length', 'TPM']
    """
    print(f"\t Processing HTSeq file: {filepath}")

    df = pd.read_csv(filepath, sep="\t", header=None, names=["gene_id", "counts"])
    df = df[~df["gene_id"].str.startswith("__")]
    df['gene_id'] = df['gene_id'].apply(strip_version)
    df['counts'] = pd.to_numeric(df['counts'], errors='coerce')

    # Merge with gene lengths
    merged = pd.merge(df, gene_lengths_df, on='gene_id', how='inner')
    # TPM calculation
    merged["TPM"] = (merged["counts"] / merged["gene_length"]) / (merged["counts"] / merged["gene_length"]).sum() * 1e6
    return merged[["gene_id", "TPM"]]

def process_salmon(filepath):
    """Process Salmon output file."""
    print(f"\t Processing Salmon file: {filepath}")
    df = pd.read_csv(filepath, sep="\t")
    df = df[["Name", "TPM"]]
    df.columns = ["gene_id", "TPM"]
    df['gene_id'] = df['gene_id'].apply(strip_version)
    return df

def process_kallisto(filepath):
    """Process Kallisto output file and compute gene-level TPM by summing transcript-level TPMs."""
    print(f"\t Processing Kallisto file: {filepath}")
    df = pd.read_csv(filepath, sep="\t")

    # gene ID (ENSG) extraction from the target_id column
    df['gene_id'] = df['target_id'].apply(lambda x: x.split('|')[1] if '|' in x and len(x.split('|')) > 1 else x)

    # Sum TPM values per gene
    gene_data = df.groupby('gene_id')['tpm'].sum().reset_index()
    gene_data.rename(columns={'tpm': 'TPM'}, inplace=True)
    gene_data['gene_id'] = gene_data['gene_id'].apply(strip_version)
    return gene_data

def process_bilatticeCount(filepath):
    """Process custom tool output file."""
    print(f'\t Processing bilatticeCount file: {filepath}')
    df = pd.read_csv(filepath, sep="\t")
    df['gene_id'] = df['gene_id'].apply(strip_version)
    return df[['gene_id', 'TPM']]

def process_rt_pcr(outputs_path):
    # Load the SEQC_qPCR file
    file_path = os.path.join(outputs_path, 'GSE56457_SEQC_qPCR_GEOSub.txt')
    data = pd.read_csv(file_path, sep="\t", comment="#", skip_blank_lines=True, engine="python")

    data["RT_A_relative"] = 2 ** (-data["SEQC_RTPCR_A"])
    data["RT_B_relative"] = 2 ** (-data["SEQC_RTPCR_B"])

    data.loc[data["RT_A_relative"] == 1.0, "RT_A_relative"] = 0.0
    data.loc[data["RT_B_relative"] == 1.0, "RT_B_relative"] = 0.0

    # Keep only gene_id and transformed values
    selected_columns = data[['ensembl_gene', 'RT_A_relative', 'RT_B_relative']]
    selected_columns.columns = ['gene_id', 'RT_A_relative', 'RT_B_relative']

    return selected_columns

def process_gt_beers(outputs_path):
    """
    Processes a ground truth BEERS2 parquet file.
    
    Steps:
    - Load file
    - Filter to rows where 'run' is in `all_bias`
    - For samples 1 through 8:
        - Filter by sample
        - Group by gene and sum TPMs
    - Return a list of 8 gene-level TPM DataFrames (one for each sample)
    """
    # Load the data
    file_path = os.path.join(outputs_path, 'beers.true_TPM.parquet')
    df = pd.read_parquet(file_path)

    # Filter only those values that are in all_bias
    df_filtered = df[df['run'] ==  'all_bias']

    # List to hold resulting gene-level DataFrames
    gene_tpm_dfs = []

    # Process each sample
    for sample_num in range(1, 9):
        # Filter for the sample
        sample_df = df_filtered[df_filtered['sample'] == sample_num]

        # Group by gene and sum TPM
        gene_tpm = sample_df.groupby('GeneID')['TPM'].sum().reset_index()
        gene_tpm.rename(columns={'GeneID': 'gene_id'}, inplace=True)

        gene_tpm_dfs.append(gene_tpm)

    return gene_tpm_dfs

def process_simple_directory(directory_path, tool_type, gene_lengths_df=None):
    all_dataframes = []
    
    for filename in os.listdir(directory_path):
        if filename.endswith(".txt"):
            filepath = os.path.join(directory_path, filename)
            if os.path.getsize(filepath) == 0:
                continue
            df = pd.DataFrame()
            if tool_type == 'HTSeq':
                df = process_htseq(filepath, gene_lengths_df=gene_lengths_df)
            elif tool_type == 'featureCounts':
                df = process_featureCounts(filepath)
            elif tool_type == 'bilatticeCount':
                df = process_bilatticeCount(filepath)
            
            if 'bowtie2' in filepath:
                df = df.rename(columns={"TPM": f'bowtie2_{tool_type}_{filename[:-4]}'})  # Rename 'count' column to the filename
            elif 'star' in filepath:
                df = df.rename(columns={"TPM": f'star_{tool_type}_{filename[:-4]}'})  # Rename 'count' column to the filename

            all_dataframes.append(df)
    if all_dataframes:
        result_df = all_dataframes[0]
        for df in all_dataframes[1:]:
            result_df = pd.merge(result_df, df, on="gene_id", how="outer")
        return result_df
    else:
        return pd.DataFrame(columns=['gene_id'])  # Return an empty DataFrame if no files were processed
    
def process_complex_directory(directory_path, tool_type):
    all_dataframes = []
    for subdir in os.listdir(directory_path):
        subdir_path = os.path.join(directory_path, subdir)
        if os.path.isdir(subdir_path):
            df = pd.DataFrame()
            if tool_type == 'salmon':
                filepath = os.path.join(subdir_path, "quant.sf")
                if os.path.exists(filepath):
                    df = process_salmon(filepath)
            elif tool_type == 'kallisto':
                filepath = os.path.join(subdir_path, "abundance.tsv")
                if os.path.exists(filepath):
                    df = process_kallisto(filepath)
            if not df.empty:
                df = df.rename(columns={"TPM": f'{tool_type}_{subdir}'})
                all_dataframes.append(df)
    
    if all_dataframes:
        result_df = all_dataframes[0]
        for df in all_dataframes[1:]:
            result_df = pd.merge(result_df, df, on="gene_id", how="outer")
        return result_df
    else:
        return pd.DataFrame(columns=['gene_id'])  # Return an empty DataFrame if no files were processed

def process_seqcA(counts_path, outputs_path, gene_lengths_df=None):
    df = process_rt_pcr(outputs_path)
    df = df[['gene_id', 'RT_A_relative']]
    for dir_name in os.listdir(counts_path):
        print(f'Processing {dir_name} directory.')
        if dir_name.endswith("seqcA"):
            dir_path = os.path.join(counts_path, dir_name)
            if 'featureCounts' in dir_name:
                df_dir = process_simple_directory(dir_path, 'featureCounts')
            elif 'HTSeq' in dir_name:
                df_dir = process_simple_directory(dir_path, 'HTSeq', gene_lengths_df=gene_lengths_df)
            elif 'bilatticeCount' in dir_name:
                df_dir = process_simple_directory(dir_path, 'bilatticeCount')
            elif 'kallisto' in dir_name:
                df_dir = process_complex_directory(dir_path, 'kallisto')
            elif 'salmon' in dir_name:
                df_dir = process_complex_directory(dir_path, 'salmon')

            # Merging the processed data with the main dataframe
            df = pd.merge(df, df_dir, on="gene_id", how="left")
    df.to_csv(os.path.join(outputs_path, 'seqcA.csv'), index=False)

def process_seqcB(counts_path, outputs_path, gene_lengths_df=None):
    df = process_rt_pcr(outputs_path)
    df = df[['gene_id', 'RT_B_relative']]
    for root, dirs, files in os.walk(counts_path):
        for dir_name in dirs:
            if dir_name.endswith("seqcB"):
                dir_path = os.path.join(root, dir_name)
                if 'featureCounts' in dir_name:
                    df_dir = process_simple_directory(dir_path, 'featureCounts')
                elif 'HTSeq' in dir_name:
                    df_dir = process_simple_directory(dir_path, 'HTSeq', gene_lengths_df=gene_lengths_df)
                elif 'bilatticeCount' in dir_name:
                    df_dir = process_simple_directory(dir_path, 'bilatticeCount')
                elif 'kallisto' in dir_name:
                    df_dir = process_complex_directory(dir_path, 'kallisto')
                elif 'salmon' in dir_name:
                    df_dir = process_complex_directory(dir_path, 'salmon')

                # Merging the processed data with the main dataframe
                df = pd.merge(df, df_dir, on="gene_id", how="left")
    df.to_csv(os.path.join(outputs_path, 'seqcB.csv'), index=False)


def process_beers(counts_path, outputs_path, gene_lengths_df=None):
    # Step 1: Load ground truth TPMs per sample (list of 8 DataFrames)
    all_samples = process_gt_beers(outputs_path)  # returns sample1 to sample8



    # Step 3: Traverse all tools' directories and process files
    for root, dirs, _ in os.walk(counts_path):
        for dir_name in dirs:
            if not dir_name.endswith("beers_all_bias"):
                continue

            dir_path = os.path.join(root, dir_name)

            # Select the tool-specific processor
            if 'featureCounts' in dir_name:
                continue
                df_dir = process_simple_directory(dir_path, 'featureCounts')
            elif 'HTSeq' in dir_name:
                continue
                df_dir = process_simple_directory(dir_path, 'HTSeq', gene_lengths_df=gene_lengths_df)
            elif 'bilatticeCount' in dir_name:
                continue
                df_dir = process_simple_directory(dir_path, 'bilatticeCount')
            elif 'kallisto' in dir_name:
                df_dir = process_complex_directory(dir_path, 'kallisto')
            elif 'salmon' in dir_name:
                df_dir = process_complex_directory(dir_path, 'salmon')
            else:
                continue

            print(type(df_dir))
            print(df_dir.head())
            for col in df_dir.columns:
                if col == "gene_id":
                    continue

                match = re.search(r"S(\d+)", col)
                if not match:
                    print(f"⚠️ Skipping column (no sample number): {col}")
                    continue

                sample_number = int(match.group(1))
                if not (1 <= sample_number <= 8):
                    print(f"⚠️ Invalid sample number in column: {col}")
                    continue

                # Prepare just gene_id + this column
                temp_df = df_dir[["gene_id", col]].copy()

                # Merge into the main sample-wide table
                all_samples[sample_number - 1] = pd.merge(
                    all_samples[sample_number - 1],
                    temp_df,
                    on="gene_id",
                    how="left"
                )
                print(all_samples[sample_number - 1].head())

    # Save each sample-wide DataFrame
    for i, df in enumerate(all_samples):
        output_path = os.path.join(outputs_path, f"beers{i + 1}.csv")
        df.to_csv(output_path, sep="\t", index=False)
        print(f"✅ Saved: {output_path}")


def process_rat_directory(directory_path, tool_type, gene_lengths_df=None, main_table=None):
    all_dataframes = []
    all_mirna_dataframes = []
    
    for filename in os.listdir(directory_path):
        if filename.endswith(".txt"):
            filepath = os.path.join(directory_path, filename)
            if os.path.getsize(filepath) == 0:
                continue
            df = pd.DataFrame()
            if tool_type == 'HTSeq':
                df = process_htseq(filepath, gene_lengths_df=gene_lengths_df)
            elif tool_type == 'featureCounts':
                df = process_featureCounts(filepath)
            elif tool_type == 'bilatticeCount':
                df = process_bilatticeCount(filepath)
            
            if 'bowtie2' in filepath:
                df = df.rename(columns={"TPM": f'bowtie2_{tool_type}_{filename[:-4]}'})  # Rename 'count' column to the filename
            elif 'star' in filepath:
                df = df.rename(columns={"TPM": f'star_{tool_type}_{filename[:-4]}'})  # Rename 'count' column to the filename

            if 'mirna' in filename:
                df = pd.merge(main_table, df, left_on='mirna_ensamble_id', right_on='gene_id', how='left')
                all_mirna_dataframes.append(df)
            else:
                df = pd.merge(main_table, df, left_on='gene_ensamble_id', right_on='gene_id', how='left')
                all_dataframes.append(df)
    if all_dataframes:
        result_df = all_dataframes[0]
        result_mirna_df = all_mirna_dataframes[0]
        for df in all_dataframes[1:]:
            result_df = pd.merge(result_df, df, on="gene_id", how="outer")
        for df in all_mirna_dataframes[1:]:
            result_mirna_df = pd.merge(result_mirna_df, df, on="gene_id", how="outer")
        return result_df, result_mirna_df
    else:
        return pd.DataFrame(columns=['gene_id'])  # Return an empty DataFrame if no files were processed
    

def process_rat(counts_path, outputs_path, gene_lengths_df=None):
    main_table = pd.read_csv(os.path.join(outputs_path, 'mirna_rna_pairs.csv'), sep="\t")

    for root, dirs, _ in os.walk(counts_path):
        for dir_name in dirs:
            if not dir_name.endswith("rat"):
                continue

            dir_path = os.path.join(root, dir_name)
            df_dir = pd.DataFrame()
            # Select the tool-specific processor
            if 'featureCounts' in dir_name:
                df_dir = process_simple_directory(dir_path, 'featureCounts')
            elif 'HTSeq' in dir_name:
                continue
                #df_dir = process_simple_directory(dir_path, 'HTSeq', gene_lengths_df=gene_lengths_df)
            elif 'bilatticeCount' in dir_name:
                df_dir = process_simple_directory(dir_path, 'bilatticeCount')
            else:
                continue

            print(type(df_dir))
            print(df_dir.head())
            for col in df_dir.columns:
                if col == "gene_id":
                    continue

                match = re.search(r"S(\d+)", col)
                if not match:
                    print(f"⚠️ Skipping column (no sample number): {col}")
                    continue

                sample_number = int(match.group(1))
                if not (1 <= sample_number <= 8):
                    print(f"⚠️ Invalid sample number in column: {col}")
                    continue

                # Prepare just gene_id + this column
                temp_df = df_dir[["gene_id", col]].copy()

                # Merge into the main sample-wide table
                all_samples[sample_number - 1] = pd.merge(
                    all_samples[sample_number - 1],
                    temp_df,
                    on="gene_id",
                    how="left"
                )

    # Save each sample-wide DataFrame
    for i, df in enumerate(all_samples):
        output_path = os.path.join(outputs_path, f"rat{i + 1}.csv")
        df.to_csv(output_path, sep="\t", index=False)


if __name__ == "__main__":
    counts_path = 'genetic_data/counts'
    outputs_path = 'genetic_data/outputs'
    mus_musculus_gene_lenghts = get_gene_lengths_from_gtf('genetic_data/annotations/Mus_musculus.GRCm38.102.gtf')
    homo_sapiens_gene_lenghts = get_gene_lengths_from_gtf('genetic_data/annotations/gencode.v19.annotation.gtf')
    #process_seqcA(counts_path, outputs_path, gene_lengths_df=homo_sapiens_gene_lenghts)
    #process_seqcB(counts_path, outputs_path, gene_lengths_df=homo_sapiens_gene_lenghts)
    process_beers(counts_path, outputs_path, gene_lengths_df=mus_musculus_gene_lenghts)
