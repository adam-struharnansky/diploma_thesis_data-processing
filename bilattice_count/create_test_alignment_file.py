import pysam

def subset_bam(input_bam_path, output_bam_path, max_reads=100000):
    with pysam.AlignmentFile(input_bam_path, "rb") as infile:
        with pysam.AlignmentFile(output_bam_path, "wb", header=infile.header) as outfile:
            for i, read in enumerate(infile):
                if i >= max_reads:
                    break
                outfile.write(read)

input_bam = 'genetic_data/alignments/bowtie2_seqcA_all/SRR896663_sorted_by_name.bam'
output_bam = "genetic_data/bilattice_count/test_samples/filtered_alignments.bam"
subset_bam(input_bam, output_bam)
