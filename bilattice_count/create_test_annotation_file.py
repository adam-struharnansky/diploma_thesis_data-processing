import pandas as pd

all_gene_ids = set()

# Get all gene IDs
for chunk in pd.read_csv(
        'genetic_data/annotations/gencode.v19.annotation.gtf', 
        sep='\t', 
        comment='#', 
        header=None,
        usecols=[2, 8],  # Only need feature type and attributes for this step
        dtype={2: str, 8: str},
        chunksize=50000
    ):
    gene_rows = chunk[chunk[2] == 'gene']
    for attr_string in gene_rows[8]:
        for entry in attr_string.split(';'):
            entry = entry.strip()
            if entry.startswith('gene_id'):
                gene_id = entry.split(' ')[1].replace('"', '')
                all_gene_ids.add(gene_id)
                break

# Deterministically select the first 100 gene_ids (sorted alphabetically)
selected_gene_ids = set(sorted(all_gene_ids)[:100])

# Function to extract gene_id from attributes column
def extract_gene_id(attr_string):
    for entry in attr_string.split(';'):
        entry = entry.strip()
        if entry.startswith('gene_id'):
            return entry.split(' ')[1].replace('"', '')
    return None

# Read full file with all columns and filter
filtered_rows = []

for chunk in pd.read_csv(
        'genetic_data/annotations/gencode.v19.annotation.gtf', 
        sep='\t',
        comment='#',
        header=None,
        dtype={0: str, 1: str, 2: str, 3: int, 4: int, 5: str, 6: str, 7: str, 8: str},
        chunksize=50000
    ):
    chunk['gene_id'] = chunk[8].map(extract_gene_id)
    filtered_chunk = chunk[chunk['gene_id'].isin(selected_gene_ids)]
    filtered_rows.append(filtered_chunk.drop(columns='gene_id'))  # remove tmp column

# Concatenate and save the filtered DataFrame
result_df = pd.concat(filtered_rows, ignore_index=True)
print(result_df.shape)
result_df.to_csv('genetic_data/bilattice_count/test_samples/filtered_genes.gtf', sep='\t', header=False, index=False)
