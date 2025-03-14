from pysam import VariantFile
import numpy as np
import pandas as pd

# File paths
compressed_vcf_filename = "../Data/raw_data/ALL.chr17.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
panel_filename = "../Data/raw_data/Users/alexandrebatista/Desktop/VS_code/First_project/phase1_files/phase1_integrated_calls.20101123.ALL.panel"
exon_positions = "../Data/raw_data/Users/alexandrebatista/Desktop/VS_code/First_project/phase1_files/CARD14_exon_positions.csv"

# Lists
genotypes = []
samples = []
positions = []

# CARD14 region
chrom = "17"
start_pos = 78143829
end_pos = 78183130

# Open VCF file
with VariantFile(compressed_vcf_filename) as vcf_reader:
    try:
        for record in vcf_reader.fetch(chrom, start=start_pos, end=end_pos):
            alleles = [record.samples[x].allele_indices for x in record.samples]
            genotypes.append(alleles)
            if not samples:
                samples = list(record.samples)
            positions.append(record.pos)
    except Exception as e:
        print(f"Error fetching variants: {e}")

# Convert to NumPy array
genotypes_array = np.array(genotypes, dtype=object)

# Check array dimensions before processing
if genotypes_array.ndim == 3:
    matrix = np.count_nonzero(genotypes_array, axis=2)
else:
    print("Unexpected genotype array shape:", genotypes_array.shape)
    matrix = np.zeros((len(positions), len(samples)))

# Create DataFrame
df = pd.DataFrame(matrix.T, columns=positions, index=samples)

# Load population panel
panel_df = pd.read_csv(panel_filename, sep="\t", usecols=[0, 1, 2], 
                       names=["Samples", "Population", "Super_population"], header=None).set_index("Samples")

# Merge with population data
final_df = df.join(panel_df, how="left")

# Debug function
def debug():
    print("Final columns:", final_df.columns[-2:])
    print("Position range:", min(positions), "-", max(positions))
    print("Matrix shape:", matrix.shape)

# Categorize population
def categorize_population(df):
    df["Group"] = np.where(df["Super_population"].isin(["EUR", "AMR"]), "High", "Low")
    return df

# Apply categorization
final_df = categorize_population(final_df)

# Variant frequency calculation
def calculate_variant_frequencies(df):
    variant_frequency_dict = {"Sample": [], "Allele_frequency": [],"Super_population": [] ,"Group": []}
    variant_columns = df.columns[:-3]  # Exclude Population, Super_population, Group
    total_positions = len(variant_columns) # Total allele pairs
    total_alleles_per_sample = 2 * total_positions
    
    # Calculates the alternate allele frequency per semple and appends all info to dictionary
    for sample in df.index:
        total_alt = df.loc[sample, variant_columns].sum()
        freq = total_alt / total_alleles_per_sample
        variant_frequency_dict["Sample"].append(sample)
        variant_frequency_dict["Allele_frequency"].append(f"{freq:.6f}")
        variant_frequency_dict["Super_population"].append(df.loc[sample, "Super_population"])
        variant_frequency_dict["Group"].append(df.loc[sample, "Group"])
    
    return pd.DataFrame(variant_frequency_dict)

# Compute variant frequencies
variant_frequency_df = calculate_variant_frequencies(final_df)
