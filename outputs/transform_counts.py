import os
import pandas as pd


def strip_version(ensembl_id):
    return ensembl_id.split('.')[0]

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

def process_htseq(filepath):
    """Process HTSeq output file."""
    print(f"\t Processing HTSeq file: {filepath}")
    df = pd.read_csv(filepath, sep="\t", header=None, names=["gene_id", "tpm"])
    df = df[~df["gene_id"].str.startswith("__")]
    df.columns = ['gene_id', 'TPM']
    df['gene_id'] = df['gene_id'].apply(strip_version)
    return df

def process_salmon(filepath):
    """Process Salmon output file."""
    print(f"\t Processing Salmon file: {filepath}")
    df = pd.read_csv(filepath, sep="\t")
    df = df[["Name", "TPM"]]
    df.columns = ["gene_id", "TPM"]
    df['gene_id'] = df['gene_id'].apply(strip_version)
    return df

def process_kallisto(filepath):
    """Process Kallisto output file."""
    print(f"\t Processing Kallisto file: {filepath}")
    df = pd.read_csv(filepath, sep="\t")

    # gene ID (ENSG) extraction from the target_id column
    df['gene_id'] = df['target_id'].apply(lambda x: x.split('|')[1])

    # Aggregate raw counts and effective lengths by gene_id
    gene_counts = df.groupby('gene_id')['est_counts'].sum().reset_index()
    gene_lengths = df.groupby('gene_id')['eff_length'].sum().reset_index()

    # Merge the aggregated counts and lengths
    gene_data = pd.merge(gene_counts, gene_lengths, on='gene_id')

    # Recalculate TPM values
    gene_data['TPM'] = (gene_data['est_counts'] / gene_data['eff_length']) * 1e6 / gene_data['est_counts'].sum()
    gene_data = gene_data[['gene_id', 'TPM']]

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

def process_gt_beers():
    pass

def process_simple_directory(directory_path, tool_type):
    all_dataframes = []
    
    for filename in os.listdir(directory_path):
        if 'schweizer' in filename or 'sklar' in filename:
            continue
        if filename.endswith(".txt"):
            filepath = os.path.join(directory_path, filename)
            df = pd.DataFrame()
            if tool_type == 'HTSeq':
                df = process_htseq(filepath)
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

def process_seqcA(counts_path, outputs_path):
    df = process_rt_pcr(outputs_path)
    df = df[['gene_id', 'RT_A_relative']]
    for dir_name in os.listdir(counts_path):
        print(f'Processing {dir_name} directory.')
        if dir_name.endswith("seqcA"):
            dir_path = os.path.join(counts_path, dir_name)
            if 'featureCounts' in dir_name:
                df_dir = process_simple_directory(dir_path, 'featureCounts')
            elif 'HTSeq' in dir_name:
                df_dir = process_simple_directory(dir_path, 'HTSeq')
            elif 'bilatticeCount' in dir_name:
                df_dir = process_simple_directory(dir_path, 'bilatticeCount')
            elif 'kallisto' in dir_name:
                df_dir = process_complex_directory(dir_path, 'kallisto')
            elif 'salmon' in dir_name:
                df_dir = process_complex_directory(dir_path, 'salmon')

            # Merge the processed data with the main dataframe
            df = pd.merge(df, df_dir, on="gene_id", how="left")
    df.to_csv(os.path.join(outputs_path, 'seqcA.csv'), index=False)

def process_seqcB(counts_path, outputs_path):
    df = process_rt_pcr(outputs_path)
    df = df['gene_id', 'RT_B_relative']
    for root, dirs, files in os.walk(counts_path):
        for dir_name in dirs:
            if dir_name.endswith("seqcB"):
                dir_path = os.path.join(root, dir_name)
                if 'featureCounts' in dir_name:
                    df_dir = process_simple_directory(dir_path, 'featureCounts')
                elif 'HTSeq' in dir_name:
                    df_dir = process_simple_directory(dir_path, 'HTSeq')
                elif 'bilatticeCount' in dir_name:
                    df_dir = process_simple_directory(dir_path, 'bilatticeCount')
                elif 'kallisto' in dir_name:
                    df_dir = process_complex_directory(dir_path, 'kallisto')
                elif 'salmon' in dir_name:
                    df_dir = process_complex_directory(dir_path, 'salmon')

                # Merge the processed data with the main dataframe
                df = pd.merge(df, df_dir, on="gene_id", how="left")
    df.to_csv(os.path.join(outputs_path, 'seqcB.csv'), index=False)

def process_beers(counts_path, outputs_path):
    df = process_gt_beers(outputs_path)
    df = df['gene_id', 'beers_GT']
    for root, dirs, files in os.walk(counts_path):
        for dir_name in dirs:
            if dir_name.endswith("beers_all_bias"):
                dir_path = os.path.join(root, dir_name)
                if 'featureCounts' in dir_name:
                    df_dir = process_simple_directory(dir_path, 'featureCounts')
                elif 'HTSeq' in dir_name:
                    df_dir = process_simple_directory(dir_path, 'HTSeq')
                elif 'bilatticeCount' in dir_name:
                    df_dir = process_simple_directory(dir_path, 'bilatticeCount')
                elif 'kallisto' in dir_name:
                    df_dir = process_complex_directory(dir_path, 'kallisto')
                elif 'salmon' in dir_name:
                    df_dir = process_complex_directory(dir_path, 'salmon')

                # Merge the processed data with the main dataframe
                df = pd.merge(df, df_dir, on="gene_id", how="left")
    df.to_csv(os.path.join(outputs_path, 'beers.csv'), index=False)

if __name__ == "__main__":
    counts_path = 'genetic_data/counts'
    outputs_path = 'genetic_data/outputs'
    process_seqcA(counts_path, outputs_path)
    process_seqcB(counts_path, outputs_path)
    # process_beers(counts_path, outputs_path)
