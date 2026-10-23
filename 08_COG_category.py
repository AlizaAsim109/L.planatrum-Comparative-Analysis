

import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
from google.colab import files

# ---------------------------------------------------------
# 1. Upload file
# ---------------------------------------------------------
uploaded = files.upload()
filename = list(uploaded.keys())[0]

# ---------------------------------------------------------
# 2. Load data (supports xlsx/xls/csv/tsv)
# ---------------------------------------------------------
if filename.lower().endswith((".xlsx", ".xls")):
    df = pd.read_excel(filename)
else:
    df = pd.read_csv(filename, sep=None, engine="python", comment="#")

# ---------------------------------------------------------
# 3. Auto-detect the COG category column
# ---------------------------------------------------------
candidates = [c for c in df.columns if "COG" in c.upper() and "CAT" in c.upper()]
if not candidates:
    raise ValueError(
        f"Could not find a COG category column automatically. "
        f"Available columns: {list(df.columns)}. "
        f"Please rename the relevant column to include 'COG' and 'category'."
    )
cog_col = candidates[0]
print(f"Using column: '{cog_col}'")

# ---------------------------------------------------------
# 4. Count individual COG letters (rows may have multiple, e.g. "EG")
# ---------------------------------------------------------
letters = []
for entry in df[cog_col].dropna().astype(str):
    entry = entry.strip()
    if entry in ("", "-", "nan", "None"):
        continue
    letters.extend(list(entry))

counts = Counter(letters)

# ---------------------------------------------------------
# 5. Bar chart — all categories, alphabetically ordered
# ---------------------------------------------------------
categories = sorted(counts.keys())
frequencies = [counts[c] for c in categories]

fig, ax = plt.subplots(figsize=(14, 6))
colors = plt.cm.viridis(
    [i / max(len(categories) - 1, 1) for i in range(len(categories))]
)
ax.bar(categories, frequencies, color=colors)
ax.set_xlabel("COG Category", fontsize=12)
ax.set_ylabel("Frequency", fontsize=12)
ax.set_title("COG Functional Category Distribution", fontsize=14)
ax.grid(axis="y", alpha=0.3)
ax.set_axisbelow(True)
plt.tight_layout()
plt.savefig("COG_bar_chart.png", dpi=300)
plt.show()

# ---------------------------------------------------------
# 6. Pie chart — top 10 categories + "Other"
# ---------------------------------------------------------
ranked = counts.most_common()
top10 = ranked[:10]
other_total = sum(v for _, v in ranked[10:])

pie_labels = [k for k, v in top10]
pie_values = [v for k, v in top10]
if other_total > 0:
    pie_labels.append("Other")
    pie_values.append(other_total)

fig2, ax2 = plt.subplots(figsize=(9, 9))
ax2.pie(
    pie_values,
    labels=pie_labels,
    autopct="%1.1f%%",
    startangle=90,
    wedgeprops={"edgecolor": "white"},
)
ax2.set_title("COG Categories (Top 10)", fontsize=14)
plt.tight_layout()
plt.savefig("COG_pie_chart.png", dpi=300)
plt.show()

# ---------------------------------------------------------
# 7. Download the figures
# ---------------------------------------------------------
files.download("COG_bar_chart.png")
files.download("COG_pie_chart.png")