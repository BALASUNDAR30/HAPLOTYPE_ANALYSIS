# halpotype anlaysis after generating CSV file

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load haplotype table
df = pd.read_csv("/content/haplotype_table.csv")

# Separate SNPs and metadata
hap_ids = df["ALLELE"]
frequencies = df["freq"]
snp_data = df.drop(columns=["ALLELE", "freq"])

# Encode alleles to integers for heatmap coloring
alleles = pd.unique(snp_data.values.ravel())
allele_to_int = {allele: idx for idx, allele in enumerate(alleles)}
int_data = snp_data.replace(allele_to_int)

# Create a color palette
palette = sns.color_palette("Set3", n_colors=len(alleles))

# Plot heatmap
fig, ax = plt.subplots(figsize=(25, 10))
sns.heatmap(
    int_data,
    annot=snp_data.values,       # Show original allele letters
    fmt="s",
    cmap=palette,
    cbar=False,
    linewidths=0.5,
    linecolor='gray',
    xticklabels=snp_data.columns,
    yticklabels=hap_ids
)

# Add frequencies to the right of heatmap
for i, freq in enumerate(frequencies):
    ax.text(snp_data.shape[1] + 0.5, i + 0.5, str(freq), va='center', ha='left', fontsize=10)

# Add label for freq column
ax.text(snp_data.shape[1] + 0.5, -0.5, "freq", fontsize=12)

# Formatting
plt.title("Haplotype Analysis", fontsize=16)
plt.xlabel("SNP Positions")
plt.ylabel("Haplotypes")
plt.tight_layout()

# Save the figure
plt.savefig("haplotype_matrix_annotated.png", dpi=300)
plt.show()
