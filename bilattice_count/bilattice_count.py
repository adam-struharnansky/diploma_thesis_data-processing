
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


class IntersectionMode(str, Enum):
    STRICT = 'strict'
    ANY = 'any'


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
        self.feature = args.feature
        self.grouping_feature = args.grouping_feature
        self.sample_size = args.sample_size
        self.intersection_mode = args.intersection_mode
        self.min_alignment_score = args.min_alignment_score
        self.max_alignment_score = args.max_alignment_score
        self.verbose = args.verbose
        self.seed = args.seed


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
    parser.add_argument("--annotation_file", type=str, required=True, help="Path to the annotation file (GTF format).")
    parser.add_argument("--mappings_file", type=str, required=True, help="Path to the mappings file (BAM format).")
    parser.add_argument("--output_file", type=str, required=True, help="Name of the output file.")
    parser.add_argument("--investigated_genes", type=str, default=None, help="File with investigated genes. Default: all genes.")
    parser.add_argument("--chunk_size", type=int, default=50000, help="Size of chunks for processing data.")
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
    parser.add_argument("--feature", type=str, default="exon",
                        help="Feature type to extract from GTF (e.g., gene, transcript, exon). Default: exon.")
    parser.add_argument("--grouping_feature", type=str, default="gene_id",
                        help="Feature type to extract from GTF, by which the feature results will be grouped (e.g., gene_id, transcript_id). Default: gene_id.")
    parser.add_argument('--intersection_mode', type=IntersectionMode, default=IntersectionMode.STRICT, choices=list(IntersectionMode), 
                        help='Mode for computing the intersection of intervals.')
    parser.add_argument("--sample_size", type=int, default=100000, help="Number of reads to sample for scaling factor computation.")
    parser.add_argument('--min_alignment_score', type=int, default=None, help='Minimum alignment score that is in the mappings file.')
    parser.add_argument('--max_alignment_score', type=int, default=None, help='Maximum alignment score that is in the mappings file.')
    parser.add_argument("--verbose", action='store_true', help="Print additional information during processing.")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility. Default: random.")

    args = parser.parse_args()

    # Set seed to a random value if not provided
    if args.seed is None:
        args.seed = random.randint(0, 1000000)
    # Set the random seed for reproducibility
    random.seed(args.seed)
    np.random.seed(args.seed)
    
    # Input files validation
    validate_file(args.annotation_file, expected_suffixes=[".gtf", ".gtf.gz"])
    validate_file(args.mappings_file, expected_suffixes=[".bam", ".sam", ".bam.gz"])
    
    if args.investigated_genes:
        validate_file(args.investigated_genes, expected_suffixes=[".txt", ".csv", ".tsv", ".txt.gz", ".csv.gz", ".tsv.gz"])
    
    return args


def process_annotation_file(config):
    """
    Load and filter a GTF annotation file based on feature type and grouping feature type .

    This function reads the specified GTF file in chunks, filters rows by the given
    feature type (e.g., 'exon', 'gene', etc.), and extracts a specific grouping feature type 
    (e.g., 'gene_id' or 'transcript_id') from the attributes column. The result is
    a DataFrame containing only relevant entries with columns: seqname, start, end, strand,
    and the extracted grouping feature.

    Args:
        config (Namespace): Configuration object containing the following attributes:
            - annotation_file (str): Path to the GTF file (tab-separated).
            - chunk_size (int): Number of lines to read per chunk.
            - feature (str): Feature type to filter by (e.g., 'exon').
            - grouping_feature  (str): Grouping feature name to extract from the GTF attributes column.
            - verbose (bool): If True, print progress messages during processing.

    Returns:
        pandas.DataFrame: A DataFrame containing filtered GTF entries with selected columns.

    Raises:
        ValueError: If no entries matching the feature type and grouping feature type are found.
    """
    if config.verbose:
        print(f'Processing annotation file {config.annotation_file} (feature: {config.feature}, grouping feature: {config.grouping_feature})')

    chunks = []
    for i, chunk in enumerate(pd.read_csv(
        config.annotation_file, 
        sep='\t', 
        comment='#', 
        header=None,
        usecols=[0, 2, 3, 4, 6, 8], 
        dtype={0: str, 2: str, 3: int, 4: int, 6: str, 8: str},
        chunksize=config.chunk_size
    )):
        filtered_chunk = chunk[chunk[2] == config.feature].copy()
        if filtered_chunk.empty:
            continue

        def parse_attributes(attr_string):
            result = {}
            for entry in attr_string.split(';'):
                entry = entry.strip()
                if entry:
                    parts = entry.replace('"', '').split(' ')
                    if len(parts) >= 2:
                        key, value = parts[0], parts[1]
                        result[key] = value
            return result

        attr_dicts = filtered_chunk[8].apply(parse_attributes)

        # Extract grouping_feature (e.g., gene_id) and feature_id (e.g., exon_id or transcript_id)
        filtered_chunk['grouping_feature'] = attr_dicts.apply(lambda d: d.get(config.grouping_feature))
        filtered_chunk['feature'] = attr_dicts.apply(
            lambda d: d.get(config.feature) if config.feature in d else f"{d.get(config.grouping_feature)}_{filtered_chunk.index}"
        )

        # Drop rows missing grouping_feature (critical for grouping later)
        filtered_chunk = filtered_chunk.dropna(subset=['grouping_feature'])

        # Keep only necessary columns
        filtered_chunk = filtered_chunk[[0, 3, 4, 6, 'feature', 'grouping_feature']]
        chunks.append(filtered_chunk)

        if config.verbose:
            print(f'Processed {(i) * config.chunk_size + chunk.shape[0]} lines of annotation file')

    if not chunks:
        raise ValueError(f"No entries with feature '{config.feature}' and grouping feature type '{config.grouping_feature}' found in {config.annotation_file}.")

    gtf_df = pd.concat(chunks, ignore_index=True)
    gtf_df.columns = ['seqname', 'start', 'end', 'strand', 'feature' ,'grouping_feature']
    if config.verbose:
        print(f'Using {gtf_df.shape[0]} features of type {config.feature}')
    return gtf_df


def create_interval_trees(config, gtf_df):
    """
    Create interval trees for fast querying of genomic intervals.

    This function builds an interval tree for each chromosome using the start and end
    positions from a filtered GTF DataFrame. Each interval stores the gene ID and strand
    as its data payload. The result allows efficient overlap queries on a per-chromosome basis.

    Args:
        config (Namespace): Configuration object with a `verbose` (bool) attribute to control logging.
        gtf_df (pandas.DataFrame): DataFrame with columns ['seqname', 'start', 'end', 'strand', 'gene_id'].

    Returns:
        dict[str, IntervalTree]: A dictionary mapping chromosome names to IntervalTree objects.
    """
    if config.verbose:
        print(f'Creating Interval Trees')
    chromosomes = gtf_df['seqname'].unique()
    trees = {chromosome: IntervalTree() for chromosome in chromosomes}
    for _, row in gtf_df.iterrows():
        start = row['start'] - 1  # convertion to 0-based indexing, half-open
        end = row['end']          # TODO: check if end should be inclusive
        if start < end:
            trees[row['seqname']].add(Interval(start, end, (row['feature'], row['strand'])))
    return trees


def t_norm(config, value_1, value_2):
    """
    Apply a t-norm operation to two input values based on the configuration.

    This function computes a t-norm (triangular norm) between two values, used
    in fuzzy logic to model conjunction (logical AND).

    The specific t-norm is defined by `config.t_norm`, an instance of the `TNorm` Enum class.
    See the `TNorm` definition for the full list of supported operations.

    Args:
        config (Namespace): Configuration object containing a `t_norm` attribute of type `TNorm`.
        value_1 (float): First operand in the t-norm operation. Expected to be in [0, 1].
        value_2 (float): Second operand in the t-norm operation. Expected to be in [0, 1].

    Returns:
        float: Result of applying the selected t-norm to the input values.

    Raises:
        ValueError: If the t-norm specified in `config.t_norm` is not recognized.
    """
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
    """
    Apply a t-conorm operation to two input values based on the configuration.

    This function computes a t-conorm (triangular conorm), used in fuzzy logic to model
    disjunction (logical OR). It is derived from the selected t-norm
    by the duality transformation:
        S(a, b) = 1 - T(1 - a, 1 - b)

    The specific t-norm is defined by `config.t_norm`, an instance of the `TNorm` Enum class.
    See the `TNorm` definition for the full list of supported operations.

    Args:
        config (Namespace): Configuration object containing a `t_norm` attribute of type `TNorm`.
        value_1 (float): First operand in the t-conorm operation. Expected to be in [0, 1].
        value_2 (float): Second operand in the t-conorm operation. Expected to be in [0, 1].

    Returns:
        float: Result of applying the derived t-conorm to the input values.
    """
    return 1 - t_norm(config, 1 - value_1, 1 - value_2)


def t_norm_batch(config, a, b):
    if config.t_norm == TNorm.DRASTIC:
        # Equivalent to:
        # if a == 1 -> return b
        # elif b == 1 -> return a
        # else -> 0
        result = np.zeros_like(a)
        mask_a1 = (a == 1)
        mask_b1 = (b == 1)

        result[mask_a1] = b[mask_a1]
        result[mask_b1 & ~mask_a1] = a[mask_b1 & ~mask_a1]
        # All other values remain 0
        return result

    elif config.t_norm == TNorm.PRODUCT:
        return a * b

    elif config.t_norm == TNorm.MINIMUM:
        return np.minimum(a, b)

    elif config.t_norm == TNorm.LUKASIEWICZ:
        return np.maximum(0.0, a + b - 1)

    else:
        raise ValueError(f"Unknown t-norm type: {config.t_norm}")
    

def t_conorm_batch(config, a, b):
    return 1 - t_norm_batch(config, 1 - a, 1 - b)


def truth_t_norm(config, left, right):
    """
    Apply t-norm to the 'for' values and t-conorm to the 'against' values in a bilattice.

    Args:
        config (Namespace): Configuration object containing a `t_norm` attribute (TNorm Enum).
        left (tuple[float, float]): For value.
        right (tuple[float, float]): Against value.

    Returns:
        tuple[float, float]: Element of bilattice after operation.
    """
    return t_norm(config, left[0], right[0]), t_conorm(config, left[1], right[1])


def truth_t_conorm(config, left, right):
    """
    Apply t-conorm to the 'for' values and t-norm to the 'against' values in a bilattice.

    Args:
        config (Namespace): Configuration object containing a `t_norm` attribute (TNorm Enum).
        left (tuple[float, float]): For value.
        right (tuple[float, float]): Against value.

    Returns:
        tuple[float, float]: Element of bilattice after operation.
    """
    return t_conorm(config, left[0], right[0]), t_norm(config, left[1], right[1])


def big_truth_conorm(config, values):
    """
    Apply chained truth t-conorm operations over a list of bilattice elements.

    This function performs a reduction using `truth_t_conorm` on a list of billatice
    values represented as (for, against) tuples.

    Args:
        config (Namespace): Configuration object containing a `t_norm` attribute (TNorm Enum).
        values (list[tuple[float, float]]): List of bilattice elements as (for, against) tuples.

    Returns:
        tuple[float, float]: Final result after chaining `truth_t_conorm` over the list.
                             Returns (0, 0) if the input list is empty.
    """
    if len(values) == 0:
        return 0,0
    result = values[0]
    for value in values[1:]:
        result = truth_t_conorm(config, result, value)
    return result


def scaling_factor_computation(config):
    """
    Determine the scaling factor and minimum alignment score for normalization.

    If both `min_alignment_score` and `max_alignment_score` are provided in the config,
    they are used directly to compute the scaling factor. Otherwise, this function
    samples up to `config.sample_size` alignments from the specified BAM file to estimate
    the observed minimum and maximum alignment scores (based on the 'AS' tag), and
    calculates the scaling factor accordingly.

    Args:
        config (Namespace): Configuration object containing:
            - mappings_file (str): Path to the BAM file.
            - min_alignment_score (int or None): Minimum alignment score (optional).
            - max_alignment_score (int or None): Maximum alignment score (optional).
            - sample_size (int): Number of alignments to sample from the file.
            - verbose (bool): Whether to print progress messages.

    Returns:
        tuple[float, float]: A tuple containing:
            - scaling_factor (float): Difference between max and min alignment scores.
            - lowest_value (float): The minimum alignment score used (either provided or estimated).

    Raises:
        ValueError: If no alignment scores ('AS' tag) are found during sampling.
    """
    if config.max_alignment_score != None and config.min_alignment_score != None:
        if config.verbose:
            print(f'Using given arguments for lowest: {config.min_alignment_score} and highest: {config.max_alignment_score} values in mappings file')
        scaling_factor = config.max_alignment_score - config.min_alignment_score
        return scaling_factor, lowest_value
    if config.verbose:
        print(f'Finding Lowest and Highest Values in mappings File')
    
    lowest_value = float('inf')
    highest_value = float('-inf')

    with pysam.AlignmentFile(config.mappings_file, "rb") as mapping_file:
        for i, mapping in enumerate(mapping_file):
            if i >= config.sample_size:
                break
            if mapping.has_tag('AS'):
                value = mapping.get_tag('AS')
                lowest_value = min(lowest_value, value)
                highest_value = max(highest_value, value)
    
    if lowest_value == float('inf') or highest_value == float('-inf'):
        raise ValueError("No AS values found in the sample size provided.")

    scaling_factor = highest_value - lowest_value
    return scaling_factor, lowest_value


def process_mapping_file(config, annotations, trees, lowest_value, scaling_factor):
    """
    Process a BAM file to compute 'for' and 'against' values for each annotated feature.

    This function reads alignments from the BAM file and computes support values for
    different features (e.g., genes, exons, transcripts) using fuzzy logic operators.
    It handles both single-end and paired-end reads, applies strand-specific filtering,
    and normalizes alignment scores using a provided scaling factor.

    Args:
        config (Namespace): Configuration object specifying parameters such as strandness,
            t-norm type, paired-end handling, and verbosity.
        annotations (pandas.DataFrame): DataFrame containing the identifiers of features of interest.
        trees (dict[str, IntervalTree]): Mapping from reference names (e.g., chromosomes) to interval trees
            of annotated features.
        lowest_value (float): Minimum alignment score used for normalization.
        scaling_factor (float): Difference between maximum and minimum alignment scores.

    Returns:
        tuple[dict[str, list[float]], dict[str, list[float]]]: Two dictionaries mapping feature IDs to:
            - lists of 'for' values (evidence supporting the feature),
            - lists of 'against' values (evidence contradicting the feature).
    """
    if config.verbose:
        print(f'Processing Mapping File {config.mappings_file}')

    gene_for_values = {gene_id: [] for gene_id in annotations['feature']}
    gene_against_values = {gene_id: [] for gene_id in annotations['feature']}

    def normalize_alignment_score(mapping):
        """
        Normalize the alignment score of a single read mapping.

        Args:
            mapping (pysam.AlignedSegment): A read alignment.

        Returns:
            float: Normalized alignment score in [0, 1]. Returns 0.0 if the 'AS' tag is missing.
        """
        if not mapping.has_tag('AS'):
            return 0.0
        return (mapping.get_tag('AS') - lowest_value) / scaling_factor

    def process_pair(mapping_1, mapping_2):
        """
        Process a pair of mappings and compute a combined FOR value.

        Args:
            mapping_1 (pysam.AlignedSegment): The first alignment in the pair.
            mapping_2 (pysam.AlignedSegment or None): The second alignment in the pair, or None if unmapped.

        Returns:
            tuple[tuple[pysam.AlignedSegment], float]: The mapping(s) and their combined FOR value
                computed using the selected t-norm.
        """
        if mapping_2 is None:
            return ((mapping_1), 
                    t_norm(config, normalize_alignment_score(mapping_1), lowest_value))
        return ((mapping_1,mapping_2), 
                 t_norm(config, normalize_alignment_score(mapping_1), normalize_alignment_score(mapping_2)))        
    
    def choose_genes(trees, mapping):
        """
        Identify which features overlap with the given mapping, respecting strandness.

        Args:
            trees (dict[str, IntervalTree]): Interval trees indexed by reference name.
            mapping (pysam.AlignedSegment): A single alignment.

        Returns:
            list[str]: List of overlapping feature IDs that match the strandness criteria.
        """
        if mapping.reference_name not in trees:
            return []
        if config.intersection_mode == IntersectionMode.STRICT:
            overlapping_genes = trees[mapping.reference_name][mapping.reference_start:mapping.reference_end]
        elif config.intersection_mode == IntersectionMode.ANY:
            start_overlapping_genes = trees[mapping.reference_name][mapping.reference_start:mapping.reference_start + 1]
            end_overlapping_genes = trees[mapping.reference_name][mapping.reference_end - 1:mapping.reference_end]
            overlapping_genes = list(set(start_overlapping_genes + end_overlapping_genes))
        else:
            raise ValueError(f"Unsupported intersection mode: {config.intersection_mode}")

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
        """
        Compute AGAINST values for each FOR value using t-conorm logic.

        Each AGAINST value is the result of applying a t-conorm to all other FOR values
        in the same group.

        Args:
            for_values (list[float]): List of FOR values from a read or read pair.

        Returns:
            list[float]: List of AGAINST values, one for each FOR value.
        """
        result = [0.0] * len(for_values)
        for i in range(len(for_values)):
            for j in range(len(for_values)):
                if i != j:
                    result[i] = t_conorm(config, result[i], for_values[j])
        return result

    def process_read(mappings):
        """
        Process a group of alignments belonging to a single read (or read pair).

        Determines which annotated features the mappings support, computes normalized
        FOR and AGAINST values using fuzzy logic, and accumulates them into the overall
        feature support dictionaries.

        Handles both single-end and paired-end reads according to the configuration.

        Args:
            mappings (list[pysam.AlignedSegment]): List of alignments associated with one read name.
        """
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
                value = normalize_alignment_score(mapping)
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

            if config.verbose and (i % 100000 == 0) and (i > 0):
                print(f'First {i} reads processed')
        if mapping_buffer: # process the last read
            process_read(mapping_buffer)
        if config.verbose and i > 0:
            print(f'All {i} reads processed')
    
    # Convert to numpy
    for gene_id in gene_for_values:
        gene_for_values[gene_id] = np.array(gene_for_values[gene_id], dtype=np.float32)
        gene_against_values[gene_id] = np.array(gene_against_values[gene_id], dtype=np.float32)

    return gene_for_values, gene_against_values


def compute_results(config, annotations, gene_for_values, gene_against_values):
    """
    Compute expression counts and TPM values for annotated features.

    This function applies bilattice-based fuzzy logic to compute expression support values
    from FOR and AGAINST scores for each feature. It estimates expression counts using repeated
    necessity evaluations and derives TPM (Transcripts Per Million) values for normalization.

    The output is saved to a tab-delimited file with columns: Feature ID, Count, and TPM.

    Args:
        config (Namespace): Configuration object specifying output path, verbosity,
            and t-norm logic.
        annotations (pandas.DataFrame): DataFrame with feature annotations, containing at least
            'feature', 'grouping_feature', 'start', and 'end' columns.
        gene_for_values (dict[str, list[float]]): Mapping from feature IDs to FOR values.
        gene_against_values (dict[str, list[float]]): Mapping from feature IDs to AGAINST values.

    Returns:
        None: Results are written directly to the output file specified in `config.output_file`.
    """
    if config.verbose:
        print(f'Computing Counts')

    def batch_necessity(for_values, against_values, samples):
        
        n = len(for_values)
        
        # Step 1: Create random permutations
        matrix = np.empty((samples, 2, n), dtype=np.float32)
        for s in range(samples):
            perm = np.random.permutation(n)
            matrix[s, 0, :] = for_values[perm]
            matrix[s, 1, :] = against_values[perm]

        # Step 2: Cumulative aggregation along each row
        for i in range(1, n):
            matrix[:, 0, i] = t_norm_batch(config, matrix[:, 0, i - 1], matrix[:, 0, i])
            matrix[:, 1, i] = t_conorm_batch(config, matrix[:, 1, i - 1], matrix[:, 1, i])

        return matrix

    def necessity_distribution(for_vals, against_vals, samples=50):
        if len(for_vals) == 0:
            return 0

        matrix = batch_necessity(for_vals, against_vals, samples)
        avgs = []

        for i in range(matrix.shape[2]):  # over positions
            # shape: (samples, 2) slice
            sample_column = matrix[:, :, i]
            avgs.append(big_truth_conorm(config, sample_column))

        for i, (f_val, a_val) in enumerate(avgs):
            if f_val < a_val:
                return i

        return len(avgs)

    annotations['length_kb'] = (annotations['end'] - annotations['start']) / 1000.0

    counts = []
    for i, row in annotations.iterrows():
        feature_id = row['feature']
        count = necessity_distribution(
            gene_for_values[feature_id],
            gene_against_values[feature_id]
        )
        counts.append(count)
        if config.verbose and i % 5000 == 0 and i > 0:
            print(f'Processed {i} features')
    if config.verbose:
        print(f'All {annotations.shape[0]} features processed')

    annotations['count'] = counts
    annotations['rpk'] = annotations.apply(
        lambda row: row['count'] / row['length_kb'] if row['length_kb'] > 0 else 0,
        axis=1
    )

    rpk_sum = annotations['rpk'].sum()
    annotations['tpm'] = annotations['rpk'].apply(
        lambda val: (val / rpk_sum) * 1e6 if rpk_sum > 0 else 0
    )

    grouped_df = annotations.groupby('grouping_feature').agg({
        'count': 'sum',
        'tpm': 'sum'
    }).reset_index()

    with open(config.output_file, 'w') as file:
        file.write(f"{config.grouping_feature}\tCount\tTPM\n")
        for i, row in grouped_df.iterrows():
            file.write(f"{row['grouping_feature']}\t{int(row['count'])}\t{row['tpm']:.4f}\n")

            

def main():
    config = Config(parse_arguments())
    annotations = process_annotation_file(config)
    trees = create_interval_trees(config, annotations)
    scaling_factor, lowest_value = scaling_factor_computation(config)
    gene_for_values, gene_against_values = process_mapping_file(config, annotations, trees, lowest_value, scaling_factor)
    compute_results(config, annotations, gene_for_values, gene_against_values)
    
if __name__ == "__main__":
    main()
