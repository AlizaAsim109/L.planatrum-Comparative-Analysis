CELL1
from google.colab import files
uploaded = files.upload()

CELL2
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="whitegrid", font_scale=1.2)
plt.rcParams["figure.dpi"] = 300

CELL3
# your filename
file_name = "Book1 (1).xlsx"

df = pd.read_excel(file_name)

# Rename columns safely (prevents typo issues)
df.columns = df.columns.str.strip()

required_cols = {"Genome", "Bacteriocin", "Bacteriocin_Class"}
missing = required_cols - set(df.columns)

if missing:
    raise ValueError(f"Missing columns: {missing}")

# Drop completely empty rows
df = df.dropna(subset=["Genome", "Bacteriocin", "Bacteriocin_Class"])

print("Total bacteriocin hits:", df.shape[0])
print("Total genomes:", df["Genome"].nunique())
print("Total bacteriocin classes:", df["Bacteriocin_Class"].nunique())


CELL4
genome_counts = df.groupby("Genome").size()

plt.figure(figsize=(8, 6))
sns.histplot(
    genome_counts,
    bins=30,
    kde=True
)

plt.xlabel("Number of bacteriocins per genome")
plt.ylabel("Number of genomes")
plt.title("Bacteriocin Burden Per Genome")

plt.tight_layout()
plt.savefig("Fig2_Bacteriocins_Per_Genome.png", dpi=600, bbox_inches="tight")
plt.show()