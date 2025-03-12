
import argparse
import pandas as pd
import pysam
import random
import numpy as np

from collections import defaultdict
from intervaltree import Interval, IntervalTree


def parse_arguments():
    """Parse and return command-line arguments."""
    # TODO file existance checks
    parser = argparse.ArgumentParser(description="Gene expression counter tool.")
    parser.add_argument("--annotation_file", type=str, default="../annotations/Mus_musculus.GRCm38.102.gtf", help="Path to the annotation file (GTF format).")
    parser.add_argument("--mappings_file", type=str, default="../mappings/star/all_bias/S3_sorted_by_name.bam", help="Path to the mappings file (BAM format).")
    parser.add_argument("--output_file", type=str, default="my_output_star_drastic", help="Name of the output file.")
    parser.add_argument("--investigated_genes", type=str, default=None, help="File with investigated genes. Default: all genes.")
    parser.add_argument("--chunk_size", type=int, default=1000, help="Size of chunks for processing data.")
    # TODO sample size for scaling
    return parser.parse_args()


def process_annotation_file(annotation_file_name, chunk_size):
    print(f'Processing annotation file {annotation_file_name}')
    chunks = []

    for chunk in pd.read_csv(
        annotation_file_name, 
        sep='\t', 
        comment='#', 
        header=None,
        usecols=[0, 2, 3, 4, 8], 
        dtype={0: str, 2: str, 3: int, 4: int, 8: str},
        chunksize=chunk_size 
    ):
        # only interested in genes for now TODO transcripts
        filtered_chunk = chunk[chunk[2] == 'gene'].copy()
        # as all rows contains gene, we no longer need this information
        filtered_chunk = filtered_chunk.drop(2, axis=1)
        # we want to extract only the gene_id from the atribute feature
        filtered_chunk[8] = filtered_chunk[8].apply(lambda x: x.split(';')[0].split(' ')[1].replace('"', ''))
        chunks.append(filtered_chunk)

    # TODO - create another structure for transcripts
    # also, add here all the non-genes (ncRNA, etc) to the same structure as genes
    # transcripts - either compute with them as part of the gene (first, compute the gene, then look at the transcripts)
    # transcripts - compute with each of them as possibility for the gene (in between genes computation)
    gtf_df = pd.concat(chunks, ignore_index=True)
    gtf_df.columns = ['seqname', 'start', 'end', 'gene_id']
    return gtf_df


def create_interval_trees(gtf_df):
    print(f'Creating Interval Trees')
    chromosomes = gtf_df['seqname'].unique() # TODO - check for speed, and/or find better option
    trees = {chromosome: IntervalTree() for chromosome in chromosomes}
    for _, row in gtf_df.iterrows():
        trees[row['seqname']].add(Interval(row['start'], row['end'], row['gene_id']))
    return trees


def t_norm(value_1, value_2):
    if value_1 == 1:
        return value_2
    if value_2 == 1:
        return value_1
    return 0
    # return value_1 * value_2
    # return min(value_1, value_2)
    # return max(0.0, (value_1 + value_2) - 1)


def t_conorm(value_1, value_2):
    return 1 - t_norm(1 - value_1, 1 - value_2)
    # return value_1 + value_2 - value_1 * value_2
    # return max(value_1, value_2)
    # return min(value_1 + value_2, 1.0)


def truth_t_norm(left, right):
    return t_norm(left[0], right[0]), t_conorm(left[1], right[1])


def truth_t_conorm(left, right):
    return t_conorm(left[0], right[0]), t_norm(left[1], right[1])


def big_truth_norm(values):
    if len(values) == 0:
        return 0,0
    result = values[0]
    for value in values[1:]:
        result = truth_t_norm(result, value)
    return result


def big_truth_conorm(values):
    if len(values) == 0:
        return 0
    result = values[0]
    for value in values[1:]:
        result = truth_t_conorm(result, value)
    return result


def scaling_factor_computation(mappings_file_name):
    print(f'Finding Lowest and Highest Values in Annotation File')
    lowest_value = 255
    highest_value = -256

    with pysam.AlignmentFile(mappings_file_name, "rb") as mapping_file:
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
    return scaling_factor, lowest_value  # TODO return low & high


def process_mapping_file(mappings_file_name, annotations, trees, lowest_value, scaling_factor):
    print(f'Processing Mapping File {mappings_file_name}')

    # TODO - add another input, if they are interested in only on some genes
    # probably think, if we need to compute the distributions for all genes

    gene_for_values = {gene_id: [] for gene_id in annotations['gene_id']}
    gene_against_values = {gene_id: [] for gene_id in annotations['gene_id']}

    # TODO - different types of annotation files (SAP)

    def compute_interval_value(mapping):
        if not mapping.has_tag('AS'):
            return 0.0
        return (mapping.get_tag('AS') - lowest_value) / scaling_factor

    def process_pair(mapping_1, mapping_2):
        first_value = compute_interval_value(mapping_1)
        second_value = compute_interval_value(mapping_2)
        return t_norm(first_value, second_value)
    
    def process_read(mappings):
        # TODO - similar to featureCount -O option, whether to count reads, if on one position there are multiple genes

        # TODO sort mappings, to create pairs
        # najskor sa pozriet na chromozomy
        # TODO, pozriet si, ci je to dalej vzdy takto, alebo nie
        # TODO - niekedy su aj tri, alebo neparny pocet, potom to nedava uplne zmysel

        #print(len(mappings))

        values = []
        genes = []

        for i in range(len(mappings) // 2):
            
            # prejst si vsetky prvy krat, a vypocitat si tieto hodnoty
            values.append(process_pair(mappings[i], mappings[i + len(mappings) // 2]))

            if mappings[i].is_unmapped or mappings[i].reference_name not in trees:
                if mappings[i + len(mappings) // 2].is_unmapped or mappings[i + len(mappings) // 2].reference_name not in trees:
                    gene = set()
                else:
                    gene = trees[mappings[i + len(mappings) // 2].reference_name][mappings[i + len(mappings) // 2].reference_start:mappings[i + len(mappings) // 2].reference_end]
            else:
                gene = trees[mappings[i].reference_name][mappings[i].reference_start:mappings[i].reference_end]
                
            if not gene:
                genes.append(None)
            else:
                genes.append(next(iter(gene)))

        for i in range(len(mappings) // 2):
            if not genes[i]:
                continue
            gene_for_values[genes[i].data].append(values[i])
            max_against = 0.0
            for j in range(len(mappings) // 2):
                if i != j:
                    max_against = max(max_against, values[j])
            gene_against_values[genes[i].data].append(max_against)
    i = 0
    with pysam.AlignmentFile(mappings_file_name, "rb") as mapping_file:
        previous_query_name = ''
        mapping_buffer = []
        for mapping in mapping_file:
            if mapping.reference_name not in trees:
                # TODO - think about it, and maybe add it only to the negative parts of current query_name
                continue

            if mapping.query_name != previous_query_name:
                process_read(mapping_buffer)
                mapping_buffer = [mapping]
                previous_query_name = mapping.query_name
            else:
                mapping_buffer.append(mapping)

            i += 1
            
            if i % 1000000 == 0:
                print(f'First {i} processed')
    return gene_for_values, gene_against_values      


def compute_results(output_file_name, annotations, gene_for_values, gene_against_values):
   
    print(f'Computing Counts')

    

    def necessity(bilattice_values):
        permutation = random.sample(range(len(bilattice_values)), len(bilattice_values))
        result = [bilattice_values[permutation[0]]] * len(bilattice_values)
        for i in range(1, len(bilattice_values)):
            result[i] = truth_t_norm(result[i - 1], bilattice_values[permutation[i]])
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
        for _ in range(120):
            values.append(necessity(billatice_values))
        for value_sample in range(len(billatice_values)):
            column = [row[value_sample] for row in values]
            avgs.append(big_truth_conorm(column))

        for i in range(len(avgs)):
            if avgs[i][0] < avgs[i][1]:
                return i
        return len(avgs)


    with open(output_file_name, 'w') as file:
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
    args = parse_arguments()
    annotations = process_annotation_file(args.annotation_file, args.chunk_size)
    trees = create_interval_trees(annotations)
    scaling_factor, lowest_value = scaling_factor_computation(args.mappings_file)
    gene_for_values, gene_against_values = process_mapping_file(args.mappings_file, annotations, trees, lowest_value, scaling_factor)
    compute_results(args.output_file, annotations, gene_for_values, gene_against_values)


if __name__ == "__main__":
    main()
