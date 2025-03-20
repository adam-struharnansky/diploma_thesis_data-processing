import os
import pandas as pd


def process_featurecounts(filepath):
    """Process featureCounts output file."""
    print(f"Processing featureCounts file: {filepath}")
    df = pd.read_csv(filepath, sep="\t", comment="#")
    df = df.iloc[:, [0, -1]]  # Keep gene ID and count column
    df.columns = ["gene_id", "count"]
    return df

def process_htseq(filepath):
    """Process HTSeq output file."""
    print(f"Processing HTSeq file: {filepath}")
    df = pd.read_csv(filepath, sep="\t", header=None, names=["gene_id", "count"])
    df = df[~df["gene_id"].str.startswith("__")]
    return df

def process_salmon(filepath):
    """Process Salmon output file."""
    print(f"Processing Salmon file: {filepath}")
    df = pd.read_csv(filepath, sep="\t")
    df = df[["Name", "TPM"]]
    df.columns = ["gene_id", "TPM"]
    return df

def process_kallisto(filepath):
    """Process Kallisto output file."""
    print(f"Processing Kallisto file: {filepath}")
    df = pd.read_csv(filepath, sep="\t")
    df = df[["target_id", "tpm"]]  # Keep gene ID and TPM
    df.columns = ["gene_id", "TPM"]
    return df

def process_mytool(filepath):
    """Process custom tool output file."""
    print(f"Processing my own tool's file: {filepath}")
    df = pd.read_csv(filepath, sep="\t")  # Adjust based on format
    return df

def process_file(filepath):
    """Detect the tool and process accordingly."""
    filename = os.path.basename(filepath)
    
    if "featureCounts" in filename:
        return process_featurecounts(filepath)
    elif "htseq" in filename.lower():
        return process_htseq(filepath)
    elif "salmon" in filename.lower():
        return process_salmon(filepath)
    elif "kallisto" in filename.lower():
        return process_kallisto(filepath)
    elif "mytool" in filename.lower():
        return process_mytool(filepath)
    else:
        print(f"Skipping unknown file: {filename}")
        return None

def process_directory(directory, suffixes):
    """Scan a directory for count files and process them."""
    results = {}

    for filename in os.listdir(directory):
        if any(filename.endswith(suffix) for suffix in suffixes):
            filepath = os.path.join(directory, filename)
            result = process_file(filepath)
            if result is not None:
                results[filename] = result

    return results

# Example usage
directory = "/path/to/count/files"
suffixes = [".txt", ".tsv"]  # Add other relevant suffixes
processed_data = process_directory(directory, suffixes)

# Print processed data (example)
for filename, df in processed_data.items():
    print(f"\nProcessed {filename}:")
    print(df.head())

def process_rt_pcr():
    # Load the SEQC_qPCR file
    file_path = 'GSE56457_SEQC_qPCR_GEOSub.txt'
    data = pd.read_csv(file_path, sep="\t", comment="#", skip_blank_lines=True, engine="python")

    data["RT_A_relative"] = 2 ** (-data["SEQC_RTPCR_A"])
    data["RT_B_relative"] = 2 ** (-data["SEQC_RTPCR_B"])

    # Keep only gene_id and transformed values
    selected_columns = data[['ensembl_gene', 'RT_A_relative', 'RT_B_relative']]


    # Display the extracted columns
    return selected_columns

def process_featurecounts()
if __name__ == "__main__":
    df = process_rt_pcr()



