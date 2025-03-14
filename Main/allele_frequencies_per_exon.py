import pandas as pd
import altair as alt
from scipy.stats import kstest, mannwhitneyu
from exon_data_processing import fetch_variant_data, process_genotype_data

# Fetch variant data
genotypes_dic, samples_dic, positions_dic = fetch_variant_data()
# Import DataFrames
dataframes, dataframes_pop_info = process_genotype_data(genotypes_dic, samples_dic, positions_dic)

# Create a list to store data for high and low groups
grouped_data = []

# Process each exon to calculate high and low group frequencies
for exon_key, df in dataframes.items():
    pop_info = dataframes_pop_info[exon_key]
    
    # Define high and low groups based on Super_population
    high_group_mask = pop_info['Super_population'].isin(['EUR', 'EAS'])
    low_group_mask = pop_info['Super_population'].isin(['AFR', 'AMR'])
    
    # Calculate frequency for high group
    if any(high_group_mask):
        # Subset df using only the indices
        high_df = df.loc[high_group_mask[high_group_mask].index]
        # High incidence population alternate allele frequency calculations
        high_alt_alleles = high_df.sum().sum()
        high_total_alleles = high_df.shape[0] * high_df.shape[1] * 2
        high_freq = high_alt_alleles / high_total_alleles if high_total_alleles > 0 else 0
        # Appends the high group information to grouped_data list
        grouped_data.append({"Exon": exon_key, "Group": "High", "Allele_Frequency": high_freq})
    
    # Calculate frequency for low group
    if any(low_group_mask):
        # Correctly subset df using only the indices where the mask is True
        low_df = df.loc[low_group_mask[low_group_mask].index]
        # Low incidence population alternate allele frequency calculations
        low_alt_alleles = low_df.sum().sum()
        low_total_alleles = low_df.shape[0] * low_df.shape[1] * 2
        low_freq = low_alt_alleles / low_total_alleles if low_total_alleles > 0 else 0
        # Appends the low group information to grouped_data list
        grouped_data.append({"Exon": exon_key, "Group": "Low", "Allele_Frequency": low_freq})

# Convert to DataFrame
grouped_freq_df = pd.DataFrame(grouped_data)

# Extract exon number and sort numerically
grouped_freq_df["Exon_Number"] = grouped_freq_df["Exon"].str.extract("(\d+)").astype(int)
grouped_freq_df = grouped_freq_df.sort_values(by="Exon_Number")

# Get the ordered list of exons
exon_order = sorted(grouped_freq_df["Exon"].unique(), 
                    key=lambda x: int(x.replace("Exon", "")))


# Plotting section
y_max_value = grouped_freq_df["Allele_Frequency"].max() # Find the maximum value in your dataset

# Create grouped bar chart
grouped_bar_chart = alt.Chart(grouped_freq_df).mark_bar().encode(
    x=alt.X("Exon:N", title="Exon", sort=exon_order),
    y=alt.Y("Allele_Frequency:Q", title="Mean Alternative Allele Frequency",
    scale=alt.Scale(domain=(0, y_max_value))), # y axis limits
    xOffset="Group:N",  # This creates side-by-side bars within each exon
    color=alt.Color("Group:N", 
                    scale=alt.Scale(domain=["High", "Low"], 
                                    range=["#1f77b4", "#ff7f0e"]))
).properties(
    title="Mean Alternative Allele Frequency per Exon by Population Group",
    width=600,
    height=400
)

# Calculate differences between high and low groups
diff_data = []
for exon in exon_order:
    exon_data = grouped_freq_df[grouped_freq_df["Exon"] == exon]
    if len(exon_data) == 2:  # Both high and low exist
        high_value = exon_data[exon_data["Group"] == "High"]["Allele_Frequency"].values[0]
        low_value = exon_data[exon_data["Group"] == "Low"]["Allele_Frequency"].values[0]
        diff = high_value - low_value
        percent_diff = (diff / max(low_value, 1e-6)) * 100  # Avoid division by zero
        diff_data.append({
            "Exon": exon,
            "Difference": diff,
            "Percent_Difference": percent_diff
        })

# Convert difference data to DataFrame if any differences were calculated
if diff_data:
    diff_df = pd.DataFrame(diff_data)
    
    # Create text labels for differences
    text_labels = alt.Chart(diff_df).mark_text(
        align='center',
        baseline='bottom',
        dy=-5,
        fontSize=10
    ).encode(
        x=alt.X("Exon:N", sort=exon_order),
        y=alt.Y("Difference:Q", title=""),
        color=alt.value("black")
    )
    
    # Combine the charts
    final_chart = grouped_bar_chart + text_labels  # Layer the text on top of the bars
    
    #Display Chart
    final_chart.display()
else:
    # Display the original chart if no differences
    grouped_bar_chart.display()
    
    
# Statistical analysis

# Normality test
print("Kruskal-Wallis test")
allele_frequencies = pd.DataFrame(grouped_data)["Allele_Frequency"].dropna()

# Perform KS test against normal distribution (using sample mean & std)
ks_stat, p_value = kstest(allele_frequencies, "norm", args=(allele_frequencies.mean(), allele_frequencies.std()))

print(f"KS Statistic: {ks_stat:.4f}, p-value= {p_value:.4f}")

# Interpretation
if p_value < 0.05:
    print("Data is NOT normally distributed (reject H0).")
else:
    print("Data may be normally distributed (fail to reject H0).")
    

# # Perform the Mann-Whitney U test for each unique exon
print("\nMann-Whitney test")
for exon in grouped_freq_df["Exon"].unique():
    exon_data = grouped_freq_df[grouped_freq_df["Exon"] == exon]
    
    # Extract allele frequency values for "High" and "Low" groups
    high_values = exon_data[exon_data["Group"] == "High"]["Allele_Frequency"]
    low_values = exon_data[exon_data["Group"] == "Low"]["Allele_Frequency"]
    
    # Perform the Mann-Whitney U test (non-parametric test for comparing distributions)
    u_stat, p_value = mannwhitneyu(high_values, low_values, alternative="two-sided")
    print(f"\nExon {exon}: U-stat={u_stat:.3f}, p-value={p_value:.4f}")

    # Interpretation

    if p_value < 0.05:
        print(f"  → Exon {exon}: Significant difference between High and Low groups.")
    else:
        print(f"  → Exon {exon}: No significant difference between High and Low groups.")