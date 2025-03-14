from pysam import VariantFile  # type: ignore
import numpy as np
import pandas as pd
from sklearn import decomposition

# File paths
compressed_vcf_filename = "/Users/alexandrebatista/Desktop/VS_code/First_project/Data/raw_data/ALL.chr17.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
panel_filename = "/Users/alexandrebatista/Desktop/VS_code/First_project/Data/raw_data/phase1_integrated_calls.20101123.ALL.panel"
exon_positions = "/Users/alexandrebatista/Desktop/VS_code/First_project/Data/raw_data/CARD14_exon_positions.csv"

# Load Exon positions
chrom = "17"
exon_positions_df = pd.read_csv(exon_positions, sep=";")
exons = [(row["Start_Position"], row["End_Position"]) for index, row in exon_positions_df.iloc[3:23].iterrows()]

# Function to fetch genotypes and positions
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

# Function to process genotype data into a DataFrame
def process_genotype_data(genotypes_dic, samples_dic, positions_dic):
    panel_df = pd.read_csv(panel_filename, sep="\t", usecols=[0, 1, 2], 
                           names=["Samples", "Population", "Super_population"], header=None)
    
    genotypes_arrays = {}
    matrices = {}
    new_df = {}

    for exon_key in genotypes_dic.keys():
        try:
            genotypes_arrays[exon_key] = np.array(genotypes_dic[exon_key])
            matrices[exon_key] = np.count_nonzero(genotypes_arrays[exon_key], axis=2)

            columns_with_prefix = [f"{exon_key}_{pos}" for pos in positions_dic[exon_key]]
            new_df[exon_key] = pd.DataFrame(
                matrices[exon_key].T, 
                columns=columns_with_prefix, 
                index=samples_dic[exon_key]
            )
        except Exception as e:
            print(f"Processing failed for {exon_key}: {e}")

    all_samples_df = pd.concat(new_df.values(), axis=1)

    combined_dataframes_pop = pd.merge(
        all_samples_df,
        panel_df,
        left_index=True,
        right_on="Samples",
        how="inner"
    ).set_index("Samples")

    # Categorizing into High and Low incidence groups
    combined_dataframes_pop["Group"] = np.nan
    # Ensure 'Group' column is treated as a string before assignment
    combined_dataframes_pop["Group"] = combined_dataframes_pop["Group"].astype("object")
    combined_dataframes_pop.loc[combined_dataframes_pop["Super_population"].isin(["EUR", "AMR"]), "Group"] = "High"
    combined_dataframes_pop.loc[combined_dataframes_pop["Super_population"].isin(["AFR", "ASN"]), "Group"] = "Low"

    return combined_dataframes_pop

# PCA Function
def perform_pca(df):
    
    exon_df_snps = df.iloc[:,:-3]  # Extract only SNP data (exclude metadata)
    matrix = exon_df_snps.to_numpy()

    pca = decomposition.PCA(n_components=2)
    pca.fit(matrix)

    transformed_data = pca.transform(matrix)
    
    # Returns % variance explained
    explained_variance = pca.explained_variance_ratio_  

    PC1_variance = print(f"PC1 explains {explained_variance[0] * 100:.2f}% of variance")
    PC2_variance = print(f"PC2 explains {explained_variance[1] * 100:.2f}% of variance")

    return transformed_data, PC1_variance, PC2_variance
