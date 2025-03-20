
import argparse
import numpy as np
import os
import pandas as pd
import pysam
import random
import sys
import warnings

from collections import defaultdict
from enum import Enum
from intervaltree import Interval, IntervalTree


class TNorm(str, Enum):
    DRASTIC = "drastic"
    PRODUCT = "product"
    MINIMUM = "minimum"
    LUKASIEWICZ = "Łukasiewicz"


class Strandness(str, Enum):
    UNSTRANDED = "unstranded"
    STRANDED = "stranded"
    REVERSE = "reverse"


class UnmappedMateHandling(str, Enum):
    ZERO = "zero"
    DISCARD = "discard"


class Config:
    def __init__(self, args):
        self.annotation_file = args.annotation_file
        self.mappings_file = args.mappings_file
        self.output_file = args.output_file
        self.investigated_genes = args.investigated_genes
        self.chunk_size = args.chunk_size
        self.t_norm = args.t_norm
        self.t_norm_param = args.t_norm_param
        self.strandness = args.strandness
        self.paired_end = args.paired_end
        self.handle_unmapped_mate = args.handle_unmapped_mate
        self.verbose = args.verbose


def validate_file(path, expected_suffixes):
    """
    Check if a file exists and has one of the expected suffixes.

    Args:
        path (str): Path to the file to validate.
        expected_suffixes (list[str] or tuple[str]): List or tuple of valid suffixes (e.g. ['.bam', '.gtf']).

    Raises:
        FileNotFoundError: If the file does not exist.
        UserWarning: If the file exists but does not have one of the expected suffixes.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Error: File '{path}' not found.")
    
    if not any(path.endswith(suffix) for suffix in expected_suffixes):
        warnings.warn(f"Warning: File '{path}' does not have an expected suffix {expected_suffixes}.")


def parse_arguments():
    """
    Parse and return command-line arguments for the gene expression counter tool.

    This function sets up an argument parser with various parameters related to input files,
    t-norm selection, strand-specific counting, and other runtime options. It also validates
    the existence and suffix of the specified input files.

    Returns:
        argparse.Namespace: Parsed arguments from the command line.

    Raises:
        FileNotFoundError: If any required input file does not exist.
        UserWarning: If any input file does not have an expected suffix.
    """

    parser = argparse.ArgumentParser(description="Gene expression counter tool.")
    parser.add_argument("--annotation_file", type=str, default="genetic_data/annotations/Mus_musculus.GRCm38.102.gtf", help="Path to the annotation file (GTF format).")
    parser.add_argument("--mappings_file", type=str, default="genetic_data/mappings/star/all_bias/S3_sorted_by_name.bam", help="Path to the mappings file (BAM format).")
    parser.add_argument("--output_file", type=str, default="genetic_data/outputs", help="Name of the output file.")
    parser.add_argument("--investigated_genes", type=str, default=None, help="File with investigated genes. Default: all genes.")
    parser.add_argument("--chunk_size", type=int, default=1000, help="Size of chunks for processing data.")
    parser.add_argument("--t_norm", type=TNorm, choices=list(TNorm), default=TNorm.PRODUCT,
                        help="Type of t-norm to use. Choices: drastic, product, minimum, Łukasiewicz.")
    parser.add_argument("--t_norm_param", type=float, default=None, 
                        help="Optional numerical parameter for certain t-norms (if applicable).")
    parser.add_argument("--strandness", type=Strandness, choices=list(Strandness), default=Strandness.UNSTRANDED,
                        help="Strand-specific counting mode. Choices: unstranded, stranded, reverse.")
    parser.add_argument("--paired_end", action='store_true', help="Set if reads are paired-end.")
    parser.add_argument("--handle_unmapped_mate", type=UnmappedMateHandling,
                        choices=list(UnmappedMateHandling), default=UnmappedMateHandling.DISCARD,
                        help="How to handle cases where only one mate is mapped: 'zero' to use 0 as score, 'discard' to skip.")
    parser.add_argument("--verbose", action='store_true', help="Print additional information during processing.")

    args = parser.parse_args()
    
    # Input files validation
    validate_file(args.annotation_file, expected_suffixes=[".gtf", ".gtf.gz"])
    validate_file(args.mappings_file, expected_suffixes=[".bam", ".sam", ".bam.gz"])
    
    if args.investigated_genes:
        validate_file(args.investigated_genes, expected_suffixes=[".txt", ".csv", ".tsv", ".txt.gz", ".csv.gz", ".tsv.gz"])
    
    return args


def process_annotation_file(config):
    if config.verbose:
        print(f'Processing annotation file {config.annotation_file} (feature: {config.feature_type}, attribute: {config.attribute})')

    chunks = []
    for chunk in pd.read_csv(
        config.annotation_file, 
        sep='\t', 
        comment='#', 
        header=None,
        usecols=[0, 2, 3, 4, 6, 8], 
        dtype={0: str, 2: str, 3: int, 4: int, 6: str, 8: str},
        chunksize=config.chunk_size
    ):
        # Filter rows by selected feature type
        filtered_chunk = chunk[chunk[2] == config.feature_type].copy()
        if filtered_chunk.empty:
            continue

        # Drop feature type column
        filtered_chunk = filtered_chunk.drop(2, axis=1)

        # Extract the attribute (e.g., gene_id or transcript_id)
        def extract_attribute(attr_string):
            for entry in attr_string.split(';'):
                entry = entry.strip()
                if entry.startswith(config.attribute):
                    return entry.split(' ')[1].replace('"', '')
            return None  # in case attribute is not found

        filtered_chunk[8] = filtered_chunk[8].apply(extract_attribute)
        filtered_chunk = filtered_chunk.dropna(subset=[8])
        chunks.append(filtered_chunk)

    if not chunks:
        raise ValueError(f"No entries with feature '{config.feature_type}' and attribute '{config.attribute}' found in {config.annotation_file}.")

    gtf_df = pd.concat(chunks, ignore_index=True)
    gtf_df.columns = ['seqname', 'start', 'end', 'strand', config.attribute]
    return gtf_df

def create_interval_trees(gtf_df):
    print(f'Creating Interval Trees')
    chromosomes = gtf_df['seqname'].unique()
    trees = {chromosome: IntervalTree() for chromosome in chromosomes}
    for _, row in gtf_df.iterrows():
        trees[row['seqname']].add(Interval(row['start'], row['end'], (row['gene_id'], row['strand'])))
    return trees


def t_norm(config, value_1, value_2):
    if config.t_norm == TNorm.DRASTIC:
        if value_1 == 1:
            return value_2
        if value_2 == 1:
            return value_1
        return 0
    elif config.t_norm == TNorm.PRODUCT:
        return value_1 * value_2
    elif config.t_norm == TNorm.MINIMUM:
        return min(value_1, value_2)
    elif config.t_norm == TNorm.LUKASIEWICZ:
        return max(0.0, (value_1 + value_2) - 1)
    else:
        raise ValueError(f"Unknown t-norm type: {config.t_norm}")


def t_conorm(config, value_1, value_2):
    return 1 - t_norm(config, 1 - value_1, 1 - value_2)


def truth_t_norm(config, left, right):
    return t_norm(config, left[0], right[0]), t_conorm(config, left[1], right[1])


def truth_t_conorm(config, left, right):
    return t_conorm(config, left[0], right[0]), t_norm(config, left[1], right[1])


def big_truth_norm(config, values):
    if len(values) == 0:
        return 0,0
    result = values[0]
    for value in values[1:]:
        result = truth_t_norm(config, result, value)
    return result


def big_truth_conorm(config, values):
    if len(values) == 0:
        return 0
    result = values[0]
    for value in values[1:]:
        result = truth_t_conorm(config, result, value)
    return result


def scaling_factor_computation(config):
    print(f'Finding Lowest and Highest Values in Annotation File')
    lowest_value = 255
    highest_value = -256

    with pysam.AlignmentFile(config.mappings_file, "rb") as mapping_file:
        sample_size = 100000
        for mapping in mapping_file:
            if mapping.has_tag('AS'):

                value = mapping.get_tag('AS')
                lowest_value = min(lowest_value, value)
                highest_value = max(highest_value, value)
            sample_size -= 1
            if sample_size <= 0:
                break

    scaling_factor = highest_value - lowest_value
    return scaling_factor, lowest_value


def process_mapping_file(config, annotations, trees, lowest_value, scaling_factor):
    print(f'Processing Mapping File {config.mappings_file}')

    gene_for_values = {gene_id: [] for gene_id in annotations['gene_id']}
    gene_against_values = {gene_id: [] for gene_id in annotations['gene_id']}

    def compute_interval_value(mapping):
        if not mapping.has_tag('AS'):
            return 0.0
        return (mapping.get_tag('AS') - lowest_value) / scaling_factor

    def process_pair(mapping_1, mapping_2):
        '''
        Process a pair of mappings (paired-end). It return a tuple with the mappings and computed FOR value of the pair.
        '''
        if mapping_2 is None:
            return ((mapping_1), 
                    t_norm(config, compute_interval_value(mapping_1), lowest_value))
        return ((mapping_1,mapping_2), 
                 t_norm(config, compute_interval_value(mapping_1), compute_interval_value(mapping_2)))        
    
    def choose_genes(trees, mapping):
        if mapping.reference_name not in trees:
            return []
        overlapping_genes = trees[mapping.reference_name][mapping.reference_start:mapping.reference_end]
        if not overlapping_genes:
            return []
        chosen_genes = []
        for interval in overlapping_genes:
            gene_id, strand = interval.data
            if config.strandness == Strandness.UNSTRANDED:
                chosen_genes.append(gene_id)
            elif config.strandness == Strandness.STRANDED:
                if (strand == '-' and mapping.is_reverse) or (strand == '+' and not mapping.is_reverse):
                    chosen_genes.append(gene_id)
            elif config.strandness == Strandness.REVERSE:
                if (strand == '-' and not mapping.is_reverse) or (strand == '+' and mapping.is_reverse):
                    chosen_genes.append(gene_id)
        return chosen_genes
    
    def compute_against_values(for_values):
        result = [0.0] * len(for_values)
        for i in range(len(for_values)):
            for j in range(len(for_values)):
                if i != j:
                    result[i] = t_conorm(config, result[i], for_values[j])
        return result

    def process_read(mappings):
        if config.paired_end: # paired-end logic
            mappings.sort(key=lambda x: (x.reference_name, x.reference_start))
            mapping_value_pairs = []

            i = 0
            while i < len(mappings):
                mapping = mappings[i]
                if mapping.mate_is_unmapped or i == len(mappings) - 1:                
                    if config.handle_unmapped_mate == UnmappedMateHandling.ZERO:
                        mapping_value_pairs.append(process_pair(mapping, None))
                    i += 1
                else:
                    mapping_value_pairs.append(process_pair(mapping, mappings[i + 1]))
                    i += 2

            gene_lists = []
            for_values = []
            for pair in mapping_value_pairs:
                genes = set()
                for mapping in pair[0]:
                    genes.update(choose_genes(trees, mapping))
                gene_lists.append((genes))
                for_values.append(pair[1])
            against_values = compute_against_values(for_values)

            for i in range(len(gene_lists)):
                for gene in gene_lists[i]:
                    gene_for_values[gene].append(for_values[i])
                    gene_against_values[gene].append(against_values[i])

        else: # single-end logic
            gene_lists = []
            for_values = []
            for mapping in mappings:
                if mapping.is_unmapped or mapping.reference_name not in trees:
                    continue
                value = compute_interval_value(mapping)
                genes = choose_genes(trees, mapping)
                if genes is None:
                    continue
                gene_lists.append(genes)
                for_values.append(value)
            against_values = compute_against_values(for_values)

            for i in range(len(gene_lists)):
                for gene in gene_lists[i]:
                    gene_for_values[gene].append(for_values[i])
                    gene_against_values[gene].append(against_values[i])
              
    with pysam.AlignmentFile(config.mappings_file, "rb") as mapping_file:
        previous_query_name = ''
        mapping_buffer = []
        for i, mapping in enumerate(mapping_file):
            if mapping.query_name != previous_query_name:
                process_read(mapping_buffer)
                mapping_buffer = []
                previous_query_name = mapping.query_name
            if mapping.is_unmapped:
                continue
            if config.handle_unmapped_mate == UnmappedMateHandling.DISCARD and mapping.mate_is_unmapped:
                continue
            mapping_buffer.append(mapping)

            if i % 10000 == 0: # TODO verbose
                print(f'First {i} processed')

        if mapping_buffer: # process the last read
            process_read(mapping_buffer)

    return gene_for_values, gene_against_values


def compute_results(config, annotations, gene_for_values, gene_against_values):
   
    print(f'Computing Counts')

    

    def necessity(bilattice_values):
        permutation = random.sample(range(len(bilattice_values)), len(bilattice_values))
        result = [bilattice_values[permutation[0]]] * len(bilattice_values)
        for i in range(1, len(bilattice_values)):
            result[i] = truth_t_norm(config, result[i - 1], bilattice_values[permutation[i]])
        return result

    def necessity_distribution(f, a):
        if len(f) == 0:
            return 0
        # treba to dat dokopy
        billatice_values = []
        for i in range(len(f)):
            billatice_values.append([f[i],a[i]])

        values = []
        avgs = []
        for _ in range(50):
            values.append(necessity(billatice_values))
        for value_sample in range(len(billatice_values)):
            column = [row[value_sample] for row in values]
            avgs.append(big_truth_conorm(config, column))

        for i in range(len(avgs)):
            if avgs[i][0] < avgs[i][1]:
                return i
        return len(avgs)


    with open(config.output_file, 'w') as file:
        # Write each row to the file
        # file.write(row + '\n')
        hjkl = 0
        #file.write(f'GeneID\tcount')

        for gene_name in annotations['gene_id']:
            number = necessity_distribution(gene_for_values[gene_name], gene_against_values[gene_name])
            #print(gene_name, number)
            file.write(f'{gene_name}\t{number}\n')
            hjkl += 1
            if hjkl % 5000 == 0:
                print(f'processed {hjkl} genes')

def main():
    config = Config(parse_arguments())
    annotations = process_annotation_file(config)
    trees = create_interval_trees(annotations)
    scaling_factor, lowest_value = scaling_factor_computation(config)
    gene_for_values, gene_against_values = process_mapping_file(config, annotations, trees, lowest_value, scaling_factor)
    compute_results(config, annotations, gene_for_values, gene_against_values)
    
if __name__ == "__main__":
    main()
