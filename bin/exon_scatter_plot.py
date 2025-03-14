import altair as alt
import pandas as pd
from bin.exon_scatter_plot_data_processing import fetch_variant_data, process_genotype_data, perform_pca

# Fetch variant data
genotypes_dic, samples_dic, positions_dic = fetch_variant_data()
combined_dataframes_pop = process_genotype_data(genotypes_dic, samples_dic, positions_dic)

# Perform PCA
exon_to_plot = perform_pca(combined_dataframes_pop)

# Create new DataFrame with PCA results
exon_to_plot["PC1"] = exon_to_plot[:, 0]
exon_to_plot["PC2"] = exon_to_plot[:, 1]


# print(exon_df_plot)

# Plot PCA 
alt.Chart(exon_to_plot).mark_point().encode(
    x=alt.X("PC1", title="Principal Component 1 (PC1)"),
    y=alt.Y("PC2", title="Principal Component 2 (PC2)"),
    color=alt.Color("Super_population", scale=alt.Scale(scheme="tableau10"))
).properties(
    title="PCA of Genetic Data by Super Population"
).display()