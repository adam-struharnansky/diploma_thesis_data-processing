
import gzip
import os
import pandas as pd
import re


def strip_version(ensembl_id):
    '''
    Strip versions from ensemble gene IDs (remove numbers after first .)

    Parameters:
    ensemble (str): ID of gene to be processed

    Returns:
    str: processed gene ID
    '''
    return ensembl_id.split('.')[0]


def get_transcript_gene_mapping(filepath):
    '''
    Parse a GTF file, and return a DataFrame with mappings gene_id and transcript_id

    Parameters:
    filepath (str): Path to the GTF file (.gtf or .gtf.gz)

    Returns:
    pd.DataFrame: DataFrame with columns ['transcript_id', 'gene_id']
    '''
    mapping = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] != "transcript":
                continue
            attributes = fields[8]
            gene_id = None
            transcript_id = None
            for attribute in attributes.split(";"):
                attribute = attribute.strip()
                if attribute.startswith("gene_id"):
                    gene_id = attribute.split('"')[1]
                if attribute.startswith("transcript_id"):
                    transcript_id = attribute.split('"')[1]
            if gene_id and transcript_id:
                mapping.append((transcript_id, gene_id))
    mapping_df = pd.DataFrame(mapping, columns=["transcript_id", "gene_id"])
    mapping_df['gene_id'] = mapping_df["gene_id"].apply(strip_version)
    mapping_df['transcript_id'] = mapping_df["transcript_id"].apply(strip_version)
    return mapping_df


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
    df['gene_id'] = df["gene_id"].apply(strip_version)
    return df


def process_featureCounts(filepath):
    """
    Process featureCounts output and compute TPM values.

    Parameters:
    filepath (str): Path to featureCounts count file

    Returns:
    pd.DataFrame: DataFrame with ['gene_id', 'TPM']
    """
    print(f"\t Processing featureCounts file: {filepath}")

    df = pd.read_csv(filepath, sep="\t", comment="#")
    df = df.iloc[:, [0, -2, -1]]  # ID, length, and count columns
    df.columns = ["gene_id", "length", "count"]
    df["gene_id"] = df["gene_id"].apply(strip_version)  # Version stripping
    df["TPM"] = (df["count"] / df["length"]) / (df["count"] / df["length"]).sum() * 1e6  # TPM calculation
    df = df[["gene_id", "TPM"]] 
    return df


def process_htseq(filepath, gene_lengths_df):
    """
    Process HTSeq output and compute TPM values using gene lengths.

    Parameters:
    filepath (str): Path to HTSeq count file
    gene_lengths_df (pd.DataFrame): DataFrame with 'gene_id' and 'gene_length' (in bp)

    Returns:
    pd.DataFrame: DataFrame with ['gene_id', 'TPM']
    """
    print(f"\t Processing HTSeq file: {filepath}")

    df = pd.read_csv(filepath, sep="\t", header=None, names=["gene_id", "counts"])
    df = df[~df["gene_id"].str.startswith("__")]
    df['gene_id'] = df['gene_id'].apply(strip_version)
    df['counts'] = pd.to_numeric(df['counts'], errors='coerce')
    merged = pd.merge(df, gene_lengths_df, on='gene_id', how='inner') # Merging with lengths
    merged["TPM"] = (merged["counts"] / merged["gene_length"]) / (merged["counts"] / merged["gene_length"]).sum() * 1e6 # TPM calculation
    return merged[["gene_id", "TPM"]]


def process_salmon(filepath):
    """
    Process Salmon output.

    Parameters:
    filepath (str): Path to Salmon count file

    Returns:
    pd.DataFrame: DataFrame with ['gene_id', 'TPM']
    """
    print(f"\t Processing Salmon file: {filepath}")

    df = pd.read_csv(filepath, sep="\t")
    df = df[["Name", "TPM"]]
    df.columns = ["gene_id", "TPM"]
    df['gene_id'] = df['gene_id'].apply(strip_version)
    return df


def process_kallisto(filepath):
    """
    Process Kallisto output.

    Parameters:
    filepath (str): Path to Kallisto count file

    Returns:
    pd.DataFrame: DataFrame with ['gene_id', 'TPM']
    """
    print(f"\t Processing Kallisto file: {filepath}")

    df = pd.read_csv(filepath, sep="\t")
    df['gene_id'] = df['target_id'].apply(lambda x: x.split('|')[1] if '|' in x and len(x.split('|')) > 1 else x) # transcript_id extraction
    gene_data = df.groupby('gene_id')['tpm'].sum().reset_index()
    gene_data.rename(columns={'tpm': 'TPM'}, inplace=True)
    gene_data['gene_id'] = gene_data['gene_id'].apply(strip_version)
    return gene_data


def process_bilatticeCount(filepath):
    """
    Process bilatticeCount output.

    Parameters:
    filepath (str): Path to bilatticeCount count file

    Returns:
    pd.DataFrame: DataFrame with ['gene_id', 'TPM']
    """
    print(f'\t Processing bilatticeCount file: {filepath}')

    df = pd.read_csv(filepath, sep="\t")
    df['gene_id'] = df['gene_id'].apply(strip_version)
    return df[['gene_id', 'TPM']]


def process_rt_pcr(outputs_path):
    """
    Process RT-PCR data file.

    Parameters:
    filepath (str): Path to RT-PCR file

    Returns:
    pd.DataFrame: DataFrame with ['gene_id', 'RT_A_relative', 'RT_B_relative']
    """
    # RT-PCR data file loading
    file_path = os.path.join(outputs_path, 'GSE56457_SEQC_qPCR_GEOSub.txt')
    data = pd.read_csv(file_path, sep="\t", comment="#", skip_blank_lines=True, engine="python")

    data["RT_A_relative"] = 2 ** (-data["SEQC_RTPCR_A"])
    data["RT_B_relative"] = 2 ** (-data["SEQC_RTPCR_B"])

    # setting the default values to be 0
    data.loc[data["RT_A_relative"] == 1.0, "RT_A_relative"] = 0.0
    data.loc[data["RT_B_relative"] == 1.0, "RT_B_relative"] = 0.0

    selected_columns = data[['ensembl_gene', 'RT_A_relative', 'RT_B_relative']]
    selected_columns.columns = ['gene_id', 'RT_A_relative', 'RT_B_relative']

    return selected_columns


def process_gt_beers(outputs_path):
    """
    Process ground truth from BEERS.

    Parameters:
    filepath (str): Path to file with ground truth values from BEERS.

    Returns:
    list(pd.DataFrame): list of DataFrames with ['gene_id', 'TPM']
    """
    # Data loading
    file_path = os.path.join(outputs_path, 'beers.true_TPM.parquet')
    df = pd.read_parquet(file_path)

    df_filtered = df[df['run'] ==  'all_bias']  # Keeping only all bias values
    gene_tpm_dfs = [] # List to hold resulting gene-level DataFrames

    # Processing of each sample
    for sample_num in range(1, 9):
        sample_df = df_filtered[df_filtered['sample'] == sample_num]  # Sample filtering
        gene_tpm = sample_df.groupby('GeneID')['TPM'].sum().reset_index()
        gene_tpm.rename(columns={'GeneID': 'gene_id'}, inplace=True)
        gene_tpm_dfs.append(gene_tpm)
    return gene_tpm_dfs


def process_simple_directory(directory_path, tool_type, gene_lengths_df=None):
    """
    Process simple directory, created by some quantification tool.

    Parameters:
    directory_path (str): Path to the directory.
    tool_type (str): What type of quantification tool was used. 
    gene_lenghts_df (pd.DataFrame): The dataframe, with gene_id and lenght of given gene (needed for some quantification tools)

    Returns:
    pd.DataFrame: a dataframe with results from each processed file in the directory ['gene_id', 'X'], where each file has one column X, that is named aligner_tool_filename
    """
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
    

def process_complex_directory(directory_path, tool_type, transcript_gene_mappings=None):
    """
    Process complex directory, created by some pseudoaligner tool.

    Parameters:
    directory_path (str): Path to the directory.
    tool_type (str): What type of quantification tool was used. 
    gene_lenghts_df (pd.Dataframe): The dataframe, with gene_id and lenght of given gene (needed for some quantification tools)

    Returns:
    pd.DataFrame: a dataframe with results from each processed file in the directory ['gene_id', 'X'], where each file has one column X, that is named aligner_tool_filename
    """
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
                if transcript_gene_mappings is not None and tool_type = 'salmon':
                    df = df.merge(transcript_gene_mappings, left_on="gene_id", right_on="transcript_id", how="left")
                    df["gene_id"] = df["gene_id_y"]  # use mapped gene_id
                    df = df.drop(["transcript_id", "gene_id_x", "gene_id_y"], axis=1)
                    # Aggregate if needed (many transcripts -> one gene)
                agg_columns = [col for col in df.columns if col != "gene_id"]
                df = df.groupby("gene_id")[agg_columns].sum().reset_index()
                df = df.rename(columns={"TPM": f'{tool_type}_{subdir}'})
                all_dataframes.append(df)
    if all_dataframes:
        result_df = all_dataframes[0]
        for df in all_dataframes[1:]:
            result_df = pd.merge(result_df, df, on="gene_id", how="outer")
        print(result_df.head())
        return result_df
    else:
        return pd.DataFrame(columns=['gene_id'])  # Return an empty DataFrame if no files were processed


def process_seqcA(counts_path, outputs_path, gene_lengths_df=None, transcript_gene_mappings=None):
    """
    Process data from SEQC A.
    Saves a csv file with results from each processed file in the directory ['gene_id', 'X'],
    where each file has one column X, that is named aligner_tool_filename.

    Parameters:
    counts_path (str): path to the counts directory, in which all the results from quantification tools are stored
    outputs_path (str): path to the directory, where the output will be stored
    gene_lenghts_df (pd.Dataframe): Dataframe with ['gene_id', 'lenght'] for each gene
    transcripts_gene_mappings (pd.Dataframe): Dataframe with ['gene_id', 'transcript_id'] for each gene and transcript

    Returns:
    None
    """
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
                df_dir = process_complex_directory(dir_path, 'kallisto', transcript_gene_mappings=transcript_gene_mappings)
            elif 'salmon' in dir_name:
                df_dir = process_complex_directory(dir_path, 'salmon', transcript_gene_mappings=transcript_gene_mappings)
            print(df_dir.head())
            # Merging the processed data with the main dataframe
            df = pd.merge(df, df_dir, on="gene_id", how="left")
    df.to_csv(os.path.join(outputs_path, 'seqcA.csv'), index=False)


def process_seqcB(counts_path, outputs_path, gene_lengths_df=None, transcript_gene_mappings=None):
    """
    Process data from SEQC B.
    Saves a csv file with results from each processed file in the directory ['gene_id', 'X'],
    where each file has one column X, that is named aligner_tool_filename.

    Parameters:
    counts_path (str): path to the counts directory, in which all the results from quantification tools are stored
    outputs_path (str): path to the directory, where the output will be stored
    gene_lenghts_df (pd.Dataframe): Dataframe with ['gene_id', 'lenght'] for each gene
    transcripts_gene_mappings (pd.Dataframe): Dataframe with ['gene_id', 'transcript_id'] for each gene and transcript

    Returns:
    None
    """
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
                    df_dir = process_complex_directory(dir_path, 'kallisto', transcript_gene_mappings=transcript_gene_mappings)
                elif 'salmon' in dir_name:
                    df_dir = process_complex_directory(dir_path, 'salmon', transcript_gene_mappings=transcript_gene_mappings)

                # Merging the processed data with the main dataframe
                df = pd.merge(df, df_dir, on="gene_id", how="left")
    df.to_csv(os.path.join(outputs_path, 'seqcB.csv'), index=False)


def process_beers(counts_path, outputs_path, gene_lengths_df=None, transcript_gene_mappings=None):
    """
    Process data from BEERS.
    Saves a cvs files for each sample with results from each processed file in the directory ['gene_id', 'X'],
    where each file has one column X, that is named aligner_tool_filename.

    Parameters:
    counts_path (str): path to the counts directory, in which all the results from quantification tools are stored
    outputs_path (str): path to the directory, where the output will be stored
    gene_lenghts_df (pd.Dataframe): Dataframe with ['gene_id', 'lenght'] for each gene
    transcripts_gene_mappings (pd.Dataframe): Dataframe with ['gene_id', 'transcript_id'] for each gene and transcript

    Returns:
    None
    """
    # ground truth loading
    all_samples = process_gt_beers(outputs_path) 

    # loop through all the subdirectories
    for root, dirs, _ in os.walk(counts_path):
        for dir_name in dirs:
            if not dir_name.endswith("beers_all_bias"):
                continue

            dir_path = os.path.join(root, dir_name)

            if 'featureCounts' in dir_name:
                df_dir = process_simple_directory(dir_path, 'featureCounts')
            elif 'HTSeq' in dir_name:
                df_dir = process_simple_directory(dir_path, 'HTSeq', gene_lengths_df=gene_lengths_df)
            elif 'bilatticeCount' in dir_name:
                df_dir = process_simple_directory(dir_path, 'bilatticeCount')
            elif 'kallisto' in dir_name:
                df_dir = process_complex_directory(dir_path, 'kallisto', transcript_gene_mappings=transcript_gene_mappings)
            elif 'salmon' in dir_name:
                df_dir = process_complex_directory(dir_path, 'salmon', transcript_gene_mappings=transcript_gene_mappings)
            else:
                continue

            # append proccesed dataframe to the correct result dataframe based on the sample number
            for col in df_dir.columns:
                if col == "gene_id":
                    continue

                match = re.search(r"S(\d+)", col)
                if not match:
                    print('Error, no sample number found in:', col)
                    continue

                sample_number = int(match.group(1))
                if not (1 <= sample_number <= 8):
                    print('Error, invalid number:', col)
                    continue

                temp_df = df_dir[["gene_id", col]].copy()

                # merging to the correct results table
                all_samples[sample_number - 1] = pd.merge(
                    all_samples[sample_number - 1],
                    temp_df,
                    on="gene_id",
                    how="left"
                )
                print(all_samples[sample_number - 1].head())

    # saving each results table
    for i, df in enumerate(all_samples):
        output_path = os.path.join(outputs_path, f"beers{i + 1}.csv")
        df.to_csv(output_path, sep="\t", index=False)


def process_rat_directory(directory_path, tool_type, gene_lengths_df=None, main_table=None):
    """
    Process data from SEQC, RAT.
    Saves a csv file with results from each processed file in the directory ['gene_id', 'X'],
    where each file has one column X, that is named aligner_tool_filename.

    Parameters:
    counts_path (str): path to the counts directory, in which all the results from quantification tools are stored
    outputs_path (str): path to the directory, where the output will be stored
    gene_lenghts_df (pd.Dataframe): Dataframe with ['gene_id', 'lenght'] for each gene
    transcripts_gene_mappings (pd.Dataframe): Dataframe with ['gene_id', 'transcript_id'] for each gene and transcript

    Returns:
    None
    """
    print(f'Processing {directory_path} directory')
    all_dataframes = []
    all_mirna_dataframes = []
    mirna_ids = set(main_table['mirna_ensamble_id'].dropna())
    target_ids = set(main_table['target_ensemble_id'].dropna())
    
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
                df = df.rename(columns={"TPM": f'bowtie2_{tool_type}_{filename[:-4]}'})
            elif 'star' in filepath:
                df = df.rename(columns={"TPM": f'star_{tool_type}_{filename[:-4]}'})

            if 'mirna' in filename:
                filtered_df = df[df['gene_id'].isin(mirna_ids)]
                all_mirna_dataframes.append(filtered_df)
            else:
                filtered_df = df[df['gene_id'].isin(target_ids)]
                all_dataframes.append(filtered_df)
    if all_dataframes and all_mirna_dataframes:
        result_df = all_dataframes[0]
        result_mirna_df = all_mirna_dataframes[0]
        for df in all_dataframes[1:]:
            result_df = pd.merge(result_df, df, on="gene_id", how="outer")
        for df in all_mirna_dataframes[1:]:
            result_mirna_df = pd.merge(result_mirna_df, df, on="gene_id", how="outer")
        return result_df, result_mirna_df
    else:
        return pd.DataFrame(columns=['gene_id']), pd.DataFrame(columns=['gene_id'])
    

def process_rat(counts_path, outputs_path, gene_lengths_df=None, gene_mapings=None):
    main_table = pd.read_csv(os.path.join(outputs_path, 'mirna_rna_pairs.csv'))
    mirna_output_table = main_table[['mirna_key','mirna_ensamble_id']].copy().drop_duplicates().dropna()
    rna_output_table = main_table[['target_gene','target_ensemble_id']].copy().drop_duplicates().dropna()
    rna_kallisto = main_table[['target_gene','target_ensemble_id']].copy().drop_duplicates().dropna()

    for root, dirs, _ in os.walk(counts_path):
        for dir_name in dirs:
            if not dir_name.endswith("rat"):
                continue

            dir_path = os.path.join(root, dir_name)
            df_rna = pd.DataFrame()

            if 'featureCounts' in dir_name:
                df_rna, df_mirna = process_rat_directory(dir_path, 'featureCounts', main_table=main_table)
            elif 'HTSeq' in dir_name:
                df_rna, df_mirna = process_rat_directory(dir_path, 'HTSeq', gene_lengths_df=gene_lengths_df, main_table=main_table)
            elif 'bilatticeCount' in dir_name:
                df_rna, df_mirna = process_rat_directory(dir_path, 'bilatticeCount', main_table=main_table)
            elif 'kallisto' in dir_name:
                df_rna_kallisto = process_complex_directory(dir_path, 'kallisto', transcript_gene_mappings=gene_mapings)
                rna_kallisto = pd.merge(rna_kallisto, df_rna_kallisto, left_on='target_ensemble_id', right_on='gene_id', how='left')
                continue
            else:
                continue

            mirna_output_table = pd.merge(mirna_output_table, df_mirna, left_on='mirna_ensamble_id', right_on='gene_id', how='left')
            mirna_output_table = mirna_output_table.drop(columns=['gene_id'])
            rna_output_table = pd.merge(rna_output_table, df_rna, left_on='target_ensemble_id', right_on='gene_id', how='left')
            rna_output_table = rna_output_table.drop(columns=['gene_id'])

    mirna_output_table.to_csv(os.path.join(outputs_path, 'mirna_rat.csv'), sep='\t', index=False)
    rna_output_table.to_csv(os.path.join(outputs_path, 'rna_rat.csv'), sep='\t', index=False)
    rna_kallisto.to_csv(os.path.join(outputs_path, 'rna_kallisto_rat.csv'), sep='\t', index=False)
    

if __name__ == "__main__":
    counts_path = 'genetic_data/counts'
    outputs_path = 'genetic_data/outputs'
    #mus_musculus_gene_lenghts = get_gene_lengths_from_gtf('genetic_data/annotations/Mus_musculus.GRCm38.102.gtf')
    #mus_musculus_mappings = get_transcript_gene_mapping('genetic_data/annotations/Mus_musculus.GRCm38.102.gtf')
    homo_sapiens_gene_lenghts = get_gene_lengths_from_gtf('genetic_data/annotations/gencode.v19.annotation.gtf')
    homo_sapines_mappings = get_transcript_gene_mapping('genetic_data/annotations/gencode.v19.annotation.gtf')
    #rattus_norvegicus_gene_lenghts = get_gene_lengths_from_gtf('genetic_data/annotations/Rattus_norvegicus.Rnor_5.0.77.gtf')
    #rattus_norvegicus_mappings = get_transcript_gene_mapping('genetic_data/annotations/Rattus_norvegicus.Rnor_5.0.77.gtf')
    process_seqcA(counts_path, outputs_path, gene_lengths_df=homo_sapiens_gene_lenghts, transcript_gene_mappings=homo_sapines_mappings)
    process_seqcB(counts_path, outputs_path, gene_lengths_df=homo_sapiens_gene_lenghts, transcript_gene_mappings=homo_sapines_mappings)
    #process_beers(counts_path, outputs_path, gene_lengths_df=mus_musculus_gene_lenghts, transcript_gene_mappings=mus_musculus_mappings)
    #process_rat(counts_path, outputs_path, rattus_norvegicus_gene_lenghts, rattus_norvegicus_mappings)
