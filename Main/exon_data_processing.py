from pysam import VariantFile  # type: ignore
import numpy as np
import pandas as pd


# File paths
compressed_vcf_filename = "../Data/raw_data/ALL.chr17.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
panel_filename = "../Data/raw_data/phase1_integrated_calls.20101123.ALL.panel"
exon_positions = "../Data/raw_data/CARD14_exon_positions.csv"

# Load Exon positions
chrom = "17"
exon_positions_df = pd.read_csv(exon_positions, sep=";")
exons = [(row["Start_Position"], row["End_Position"]) for index, row in exon_positions_df.iloc[3:23].iterrows()]

# Function to fetch genotypes and samples
def fetch_variant_data():
    genotypes_dic = {}  
    samples_dic = {}  
    positions_dic = {}  

    with VariantFile(compressed_vcf_filename) as vcf_reader:
        for i, (start_pos, end_pos) in enumerate(exons, start=3):
            exon_key = f"Exon{i}"
            genotypes_dic[exon_key] = []
            samples_dic[exon_key] = []
            positions_dic[exon_key] = []

            try:
                for record in vcf_reader.fetch(chrom, start=start_pos, end=end_pos):
                    genotypes_dic[exon_key].append([record.samples[x].allele_indices for x in record.samples])
                    samples_dic[exon_key] = [sample for sample in record.samples]
                    positions_dic[exon_key].append(record.pos)
            except Exception as e:
                print(f"Fetching dictionaries for {exon_key} failed due to: {e}")

    return genotypes_dic, samples_dic, positions_dic

# Function to process genotype data into DataFrames
def process_genotype_data(genotypes_dic, samples_dic, positions_dic):
    panel_df = pd.read_csv(panel_filename, sep="\t", usecols=[0, 1, 2], 
                           names=["Samples", "Population", "Super_population"], header=None)
    
    genotypes_arrays = {}
    matrices = {}
    dataframes = {}
    dataframes_pop_info = {}

    for exon_key in genotypes_dic.keys():
        try:
            genotypes_arrays[exon_key] = np.array(genotypes_dic[exon_key]) # Creates the genotypes array dict
            matrices[exon_key] = np.count_nonzero(genotypes_arrays[exon_key], axis=2) # Converts the arrays into a tuple for analysis

            dataframes[exon_key] = pd.DataFrame(
                matrices[exon_key].T, 
                columns=positions_dic[exon_key], 
                index=samples_dic[exon_key]
            )
            # Creates a neu df containing population information
            merged_df = pd.merge(
                dataframes[exon_key], 
                panel_df, 
                left_index=True, 
                right_on="Samples", 
                how="inner"
            ).set_index("Samples")

            column_order = positions_dic[exon_key] + ["Population", "Super_population"]
            dataframes_pop_info[exon_key] = merged_df[column_order]
        except Exception as e:
            print(f"Processing failed for {exon_key}: {e}")

    return dataframes, dataframes_pop_info



