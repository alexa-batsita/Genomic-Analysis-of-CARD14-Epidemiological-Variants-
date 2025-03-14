import pandas as pd
import numpy as np
import altair as alt
from Main.exon_data_processing import fetch_variant_data, process_genotype_data

# Fetch variant data
genotypes_dic, samples_dic, positions_dic = fetch_variant_data()
dataframes, dataframes_pop_info = process_genotype_data(genotypes_dic, samples_dic, positions_dic)

# Combine dataframes
combined_dataframes = pd.concat(dataframes.values(), axis=1, keys=dataframes.keys())
combined_dataframes_pop = pd.concat(dataframes_pop_info.values(), axis=1, keys=dataframes_pop_info.keys())

# Function to calculate allele frequencies
def calculate_variant_frequencies(df):
    variant_frequency_dict = {"Sample": [], "Allele_frequency": [], "Group": []}

    variant_columns = df.columns
    total_positions = len(variant_columns)
    total_alleles_per_sample = 2 * total_positions  

    for sample in df.index:
        total_alt = df.loc[sample, variant_columns].sum()
        freq = total_alt / total_alleles_per_sample if total_alleles_per_sample > 0 else np.nan
        variant_frequency_dict["Sample"].append(sample)
        variant_frequency_dict["Allele_frequency"].append(f"{freq:.6f}")
        variant_frequency_dict["Group"].append(group_series[sample])
    
    return pd.DataFrame(variant_frequency_dict)

for exon_key in dataframes.keys():
    super_population = combined_dataframes_pop[(exon_key), 'Super_population']
    group_series = pd.Series(index=super_population.index, dtype='object')
    group_series[super_population.isin(["EUR", "AMR"])] = "High"
    group_series[super_population.isin(["AFR", "ASN"])] = "Low"

# Calculate allele frequencies
variant_frequency_df = calculate_variant_frequencies(combined_dataframes)


# Convert allele frequency to float
variant_frequency_df["Allele_frequency"] = variant_frequency_df["Allele_frequency"].astype(float)

# Create Bar Chart
bar_chart = alt.Chart(variant_frequency_df).mark_bar().encode(
    x=alt.X("Group:N", title="Group"),
    y=alt.Y("mean(Allele_frequency):Q", title="Mean Alternative Allele Frequency"),
    color=alt.Color("Group", scale=alt.Scale(scheme="tableau10")),
).properties(
    title="Mean Alternative Allele Frequency Across Groups"
)

bar_chart.display()



# Save DataFrames to CSV
combined_dataframes.to_csv("exon_matrix.csv")
combined_dataframes_pop.to_csv("exons_matrix_pop.csv")
variant_frequency_df.to_csv("exons_frequencies.csv")
print("CSV files saved.")