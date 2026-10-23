CELL1
from google.colab import files
uploaded = files.upload()


CELL2
# ================================================
# L. PLANTARUM NEUROACTIVE GENE VISUALIZATION
# ================================================

# Install required packages (run in Colab)
!pip install pandas numpy matplotlib seaborn plotly scipy scikit-learn openpyxl -q

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import scipy.cluster.hierarchy as sch
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

# Set professional publication style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")
sns.set_context("paper", font_scale=1.5, rc={"lines.linewidth": 2.5})

# ================================================
# 1. LOAD AND PREPARE DATA
# ================================================

def load_and_merge_data(file_path='neuroactive (1).xlsx'):
    """
    Load and merge gene and genome metadata
    """
    # Read both sheets
    genes_df = pd.read_excel(file_path, sheet_name='genes')
    genomes_df = pd.read_excel(file_path, sheet_name='genomes')

    print(f"Genes sheet shape: {genes_df.shape}")
    print(f"Genomes sheet shape: {genomes_df.shape}")

    # Display column info
    print("\n=== Genes DataFrame Columns ===")
    print(genes_df.columns.tolist())

    print("\n=== Genomes DataFrame Columns ===")
    print(genomes_df.columns.tolist())

    # Clean column names
    genes_df.columns = [str(col).strip() for col in genes_df.columns]
    genomes_df.columns = [str(col).strip() for col in genomes_df.columns]

    return genes_df, genomes_df

# Load your data
genes_df, genomes_df = load_and_merge_data('neuroactive (1).xlsx')  # Update with your file name

# ================================================
# 2. DATA PREPROCESSING AND QUALITY CONTROL
# ================================================

def preprocess_data(genes_df, genomes_df):
    """
    Clean and prepare data for visualization
    """
    # Identify genome columns (assuming they start after Gene and Functional_tag)
    gene_info_cols = ['Gene', 'Functional_tag']
    genome_cols = [col for col in genes_df.columns if col not in gene_info_cols]

    print(f"\nFound {len(genome_cols)} genome columns")
    print(f"First 5 genome columns: {genome_cols[:5]}")

    # Create binary presence/absence matrix if needed
    # Assuming values are counts or binary (0/1)
    gene_matrix = genes_df[genome_cols].copy()

    # Convert to binary presence/absence if not already
    # gene_matrix = gene_matrix.applymap(lambda x: 1 if x > 0 else 0)

    # Add gene information back
    gene_matrix['Gene'] = genes_df['Gene']
    if 'Functional_tag' in genes_df.columns:
        gene_matrix['Functional_tag'] = genes_df['Functional_tag']

    return gene_matrix, genome_cols

gene_matrix, genome_cols = preprocess_data(genes_df, genomes_df)

# ================================================
# 3. GROUP-LEVEL ANALYSIS (Aggregating 318 genomes)
# ================================================

def create_group_level_matrix(gene_matrix, genomes_df, genome_cols):
    """
    Aggregate gene data by Source/Group to reduce dimensionality
    """
    # Create a mapping from genome column to source
    genome_to_source = {}
    for idx, row in genomes_df.iterrows():
        genome_id = str(row['Genome']).strip()
        source = str(row['Source']).strip()
        genome_to_source[genome_id] = source

    # Group by source
    sources = sorted(set(genome_to_source.values()))
    print(f"\nFound {len(sources)} unique sources: {sources}")

    # Create aggregated matrix
    aggregated_data = []

    for gene_idx, gene_row in gene_matrix.iterrows():
        gene_name = gene_row['Gene']
        gene_data = {'Gene': gene_name}

        if 'Functional_tag' in gene_row:
            gene_data['Functional_tag'] = gene_row['Functional_tag']

        for source in sources:
            # Get all genomes from this source
            source_genomes = [g for g, s in genome_to_source.items() if s == source]

            # Calculate prevalence in this source
            source_values = []
            for genome in source_genomes:
                if genome in gene_row:
                    val = gene_row[genome]
                    if pd.notna(val):
                        source_values.append(val)

            # Calculate metrics
            if source_values:
                prevalence = len([v for v in source_values if v > 0]) / len(source_values)
                mean_value = np.mean(source_values)
            else:
                prevalence = 0
                mean_value = 0

            gene_data[f'{source}_prevalence'] = prevalence
            gene_data[f'{source}_mean'] = mean_value

        aggregated_data.append(gene_data)

    aggregated_df = pd.DataFrame(aggregated_data)
    return aggregated_df, sources

aggregated_df, sources = create_group_level_matrix(gene_matrix, genomes_df, genome_cols)

# ================================================
# 4. PROFESSIONAL HEATMAP VISUALIZATION
# ================================================

def create_publication_heatmap(aggregated_df, sources, metric='prevalence'):
    """
    Create publication-ready clustered heatmap
    """
    # Prepare data for heatmap
    heatmap_data = []
    gene_names = []
    functional_tags = []

    for idx, row in aggregated_df.iterrows():
        gene_names.append(row['Gene'])
        if 'Functional_tag' in row:
            functional_tags.append(row['Functional_tag'])

        source_values = []
        for source in sources:
            col_name = f'{source}_{metric}'
            source_values.append(row[col_name])
        heatmap_data.append(source_values)

    heatmap_array = np.array(heatmap_data)

    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 14),
                           gridspec_kw={'height_ratios': [1, 0.1],
                                       'width_ratios': [1, 0.05]})

    # Main heatmap
    ax = axes[0, 0]

    # Perform hierarchical clustering
    try:
        # Cluster rows (genes)
        row_linkage = sch.linkage(heatmap_array, method='average', metric='euclidean')
        row_order = sch.dendrogram(row_linkage, no_plot=True)['leaves']

        # Cluster columns (sources)
        col_linkage = sch.linkage(heatmap_array.T, method='average', metric='euclidean')
        col_order = sch.dendrogram(col_linkage, no_plot=True)['leaves']

        # Reorder data
        clustered_data = heatmap_array[row_order, :]
        clustered_data = clustered_data[:, col_order]
        clustered_genes = [gene_names[i] for i in row_order]
        clustered_sources = [sources[i] for i in col_order]

    except:
        print("Clustering failed, using original order")
        clustered_data = heatmap_array
        clustered_genes = gene_names
        clustered_sources = sources

    # Create heatmap
    im = ax.imshow(clustered_data, aspect='auto', cmap='YlOrRd',
                   interpolation='nearest', vmin=0, vmax=1)

    # Set labels
    ax.set_xticks(np.arange(len(clustered_sources)))
    ax.set_yticks(np.arange(len(clustered_genes)))
    ax.set_xticklabels(clustered_sources, rotation=45, ha='right', fontsize=10)
    ax.set_yticklabels(clustered_genes, fontsize=8)

    ax.set_xlabel('Source/Environment', fontsize=12, fontweight='bold')
    ax.set_ylabel('Neuroactive Genes', fontsize=12, fontweight='bold')
    ax.set_title(f'L. plantarum Neuroactive Gene Prevalence by Source\n(n={len(gene_names)} genes across {len(sources)} environments)',
                fontsize=14, fontweight='bold', pad=20)

    # Add colorbar
    cbar_ax = axes[0, 1]
    cbar = plt.colorbar(im, cax=cbar_ax)
    cbar.set_label('Prevalence (0-1)', fontsize=10, fontweight='bold')

    # Add functional tag annotations if available
    if len(functional_tags) == len(gene_names):
        ax2 = axes[1, 0]
        colors = plt.cm.tab20(np.linspace(0, 1, len(set(functional_tags))))
        tag_to_color = {tag: colors[i] for i, tag in enumerate(sorted(set(functional_tags)))}

        for i, tag in enumerate(functional_tags):
            ax2.barh(i, 1, color=tag_to_color[tag], edgecolor='black')

        ax2.set_xlim(0, 1)
        ax2.set_ylim(0, len(functional_tags))
        ax2.set_yticks([])
        ax2.set_xlabel('Functional Categories', fontsize=10)
        ax2.set_title('Gene Function Legend', fontsize=10)

        # Create custom legend
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=tag_to_color[tag],
                               edgecolor='black',
                               label=tag) for tag in sorted(set(functional_tags))]
        axes[1, 1].legend(handles=legend_elements, loc='center',
                         fontsize=6, ncol=2)
        axes[1, 1].axis('off')

    plt.tight_layout()
    plt.savefig('neuroactive_gene_heatmap.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('neuroactive_gene_heatmap.png', dpi=300, bbox_inches='tight')
    plt.show()

    return clustered_data, clustered_genes, clustered_sources

# Generate heatmap
clustered_data, clustered_genes, clustered_sources = create_publication_heatmap(
    aggregated_df, sources, metric='prevalence'
)

# ================================================
# 5. INTERACTIVE HEATMAP (for exploration)
# ================================================

def create_interactive_heatmap(aggregated_df, sources):
    """
    Create interactive Plotly heatmap for detailed exploration
    """
    # Prepare data
    heatmap_values = []
    annotations = []

    for idx, row in aggregated_df.iterrows():
        gene_vals = []
        for source in sources:
            val = row[f'{source}_prevalence']
            gene_vals.append(val)

            # Create annotation text
            annotations.append(dict(
                x=source,
                y=row['Gene'],
                text=f"{val:.2f}",
                showarrow=False,
                font=dict(size=8)
            ))
        heatmap_values.append(gene_vals)

    # Create interactive heatmap
    fig = go.Figure(data=go.Heatmap(
        z=heatmap_values,
        x=sources,
        y=aggregated_df['Gene'],
        colorscale='YlOrRd',
        colorbar=dict(title="Prevalence"),
        hoverongaps=False,
        hoverinfo='x+y+z'
    ))

    # Update layout
    fig.update_layout(
        title={
            'text': "Interactive Neuroactive Gene Prevalence Heatmap",
            'y':0.95,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top',
            'font': {'size': 20}
        },
        xaxis_title="Source/Environment",
        yaxis_title="Neuroactive Genes",
        height=800,
        width=1200,
        template="plotly_white"
    )

    fig.write_html("interactive_heatmap.html")
    fig.show()

    return fig

# Generate interactive version
interactive_fig = create_interactive_heatmap(aggregated_df, sources)

# ================================================
# 6. PROFESSIONAL BAR CHARTS
# ================================================

def create_publication_barcharts(aggregated_df, sources):
    """
    Create multi-panel bar charts for key findings
    """
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Panel A: Top genes by overall prevalence
    ax1 = axes[0, 0]
    overall_prevalence = []
    for idx, row in aggregated_df.iterrows():
        gene_prev = np.mean([row[f'{source}_prevalence'] for source in sources])
        overall_prevalence.append((row['Gene'], gene_prev))

    # Sort and get top 20
    overall_prevalence.sort(key=lambda x: x[1], reverse=True)
    top_genes = [x[0] for x in overall_prevalence[:20]]
    top_values = [x[1] for x in overall_prevalence[:20]]

    bars1 = ax1.barh(range(len(top_genes)), top_values, color='steelblue')
    ax1.set_yticks(range(len(top_genes)))
    ax1.set_yticklabels(top_genes, fontsize=9)
    ax1.set_xlabel('Average Prevalence', fontsize=11, fontweight='bold')
    ax1.set_title('A) Top 20 Neuroactive Genes by Prevalence',
                 fontsize=12, fontweight='bold', pad=15)
    ax1.invert_yaxis()

    # Panel B: Source richness (number of genes per source)
    ax2 = axes[0, 1]
    source_richness = {}
    for source in sources:
        source_col = f'{source}_prevalence'
        genes_in_source = sum(aggregated_df[source_col] > 0.5)  # Threshold
        source_richness[source] = genes_in_source

    sources_sorted = sorted(source_richness.items(), key=lambda x: x[1], reverse=True)
    source_names = [x[0] for x in sources_sorted]
    richness_vals = [x[1] for x in sources_sorted]

    bars2 = ax2.bar(range(len(source_names)), richness_vals, color='coral')
    ax2.set_xticks(range(len(source_names)))
    ax2.set_xticklabels(source_names, rotation=45, ha='right', fontsize=9)
    ax2.set_ylabel('Number of Genes (prevalence > 0.5)', fontsize=11, fontweight='bold')
    ax2.set_title('B) Neuroactive Gene Richness by Source',
                 fontsize=12, fontweight='bold', pad=15)

    # Panel C: Gene prevalence distribution
    ax3 = axes[1, 0]
    all_prevalences = []
    for source in sources:
        source_col = f'{source}_prevalence'
        all_prevalences.extend(aggregated_df[source_col].values)

    ax3.hist(all_prevalences, bins=30, color='seagreen', alpha=0.7, edgecolor='black')
    ax3.set_xlabel('Gene Prevalence', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Frequency', fontsize=11, fontweight='bold')
    ax3.set_title('C) Distribution of Gene Prevalence Values',
                 fontsize=12, fontweight='bold', pad=15)
    ax3.grid(True, alpha=0.3)

    # Panel D: Source similarity dendrogram
    ax4 = axes[1, 1]
    # Calculate correlation matrix between sources
    source_corr = []
    for i, source1 in enumerate(sources):
        row = []
        for j, source2 in enumerate(sources):
            corr = np.corrcoef(aggregated_df[f'{source1}_prevalence'],
                              aggregated_df[f'{source2}_prevalence'])[0, 1]
            row.append(1 - corr)  # Convert to distance
        source_corr.append(row)

    source_corr = np.array(source_corr)

    # Perform hierarchical clustering
    linkage = sch.linkage(source_corr, method='average')
    dendro = sch.dendrogram(linkage, labels=sources, ax=ax4,
                           orientation='right', leaf_font_size=9)
    ax4.set_xlabel('Distance (1 - Correlation)', fontsize=11, fontweight='bold')
    ax4.set_title('D) Source Similarity Dendrogram',
                 fontsize=12, fontweight='bold', pad=15)

    plt.suptitle('L. plantarum Neuroactive Gene Landscape Analysis',
                fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig('neuroactive_gene_barcharts.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('neuroactive_gene_barcharts.png', dpi=300, bbox_inches='tight')
    plt.show()

# Generate bar charts
create_publication_barcharts(aggregated_df, sources)

# ================================================
# 7. FUNCTIONAL CATEGORY ANALYSIS
# ================================================

def analyze_functional_categories(genes_df, aggregated_df, sources):
    """
    Analyze genes by functional categories
    """
    if 'Functional_tag' not in genes_df.columns:
        print("No Functional_tag column found")
        return

    # Merge functional tags
    func_analysis = aggregated_df[['Gene']].copy()
    func_analysis['Functional_tag'] = genes_df['Functional_tag'].values

    # Calculate prevalence by functional category
    func_prevalence = {}
    for func_tag in func_analysis['Functional_tag'].unique():
        if pd.isna(func_tag):
            continue

        func_genes = func_analysis[func_analysis['Functional_tag'] == func_tag]['Gene'].tolist()
        func_mask = aggregated_df['Gene'].isin(func_genes)

        if sum(func_mask) == 0:
            continue

        tag_prevalence = {}
        for source in sources:
            source_col = f'{source}_prevalence'
            tag_prevalence[source] = aggregated_df.loc[func_mask, source_col].mean()

        func_prevalence[func_tag] = tag_prevalence

    # Create functional category heatmap
    func_df = pd.DataFrame(func_prevalence).T

    plt.figure(figsize=(12, 8))
    sns.heatmap(func_df, cmap='viridis', annot=True, fmt='.2f',
                linewidths=0.5, cbar_kws={'label': 'Average Prevalence'})
    plt.title('Neuroactive Functional Categories by Source',
             fontsize=14, fontweight='bold', pad=20)
    plt.xlabel('Source/Environment', fontsize=12)
    plt.ylabel('Functional Category', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig('functional_categories_heatmap.pdf', dpi=300, bbox_inches='tight')
    plt.show()

    return func_df

# Run functional analysis if tags are available
func_df = analyze_functional_categories(genes_df, aggregated_df, sources)

# ================================================
# 8. STATISTICAL SUMMARIES FOR PUBLICATION
# ================================================

def generate_statistical_summary(aggregated_df, sources):
    """
    Generate key statistics for manuscript
    """
    print("\n" + "="*60)
    print("STATISTICAL SUMMARY FOR MANUSCRIPT")
    print("="*60)

    # Overall statistics
    total_genes = len(aggregated_df)
    print(f"\n1. GENE-LEVEL ANALYSIS:")
    print(f"   - Total neuroactive genes analyzed: {total_genes}")

    # Prevalence distribution
    all_prevalences = []
    for source in sources:
        all_prevalences.extend(aggregated_df[f'{source}_prevalence'].values)

    print(f"   - Mean prevalence across all genes/sources: {np.mean(all_prevalences):.3f}")
    print(f"   - Prevalence standard deviation: {np.std(all_prevalences):.3f}")

    # Core genes (present in >90% of sources)
    core_genes = []
    for idx, row in aggregated_df.iterrows():
        high_prevalence_count = sum([row[f'{source}_prevalence'] > 0.9 for source in sources])
        if high_prevalence_count / len(sources) > 0.9:
            core_genes.append(row['Gene'])

    print(f"   - Core neuroactive genes (>90% prevalence): {len(core_genes)}")

    # Source-specific analysis
    print(f"\n2. SOURCE-LEVEL ANALYSIS:")
    for source in sources:
        source_genes = aggregated_df[aggregated_df[f'{source}_prevalence'] > 0.5]
        unique_genes = aggregated_df[aggregated_df[f'{source}_prevalence'] > 0.8]

        print(f"   - {source}:")
        print(f"     • Genes with >50% prevalence: {len(source_genes)}")
        print(f"     • Genes with >80% prevalence: {len(unique_genes)}")
        if len(unique_genes) > 0:
            print(f"     • Potential biomarkers: {', '.join(unique_genes['Gene'].head(3).tolist())}")

    # Correlation between sources
    print(f"\n3. SOURCE CORRELATIONS:")
    corr_matrix = []
    for i, source1 in enumerate(sources):
        for j, source2 in enumerate(sources[i+1:], i+1):
            corr = np.corrcoef(aggregated_df[f'{source1}_prevalence'],
                              aggregated_df[f'{source2}_prevalence'])[0, 1]
            if corr > 0.7:
                print(f"   - High correlation ({corr:.3f}) between {source1} and {source2}")

    # Save summary to CSV
    summary_df = pd.DataFrame({
        'Statistic': ['Total_Genes', 'Mean_Prevalence', 'Core_Genes'],
        'Value': [total_genes, np.mean(all_prevalences), len(core_genes)]
    })
    summary_df.to_csv('statistical_summary.csv', index=False)

    print(f"\nSummary saved to 'statistical_summary.csv'")

# Generate statistical summary
generate_statistical_summary(aggregated_df, sources)

# ================================================
# 9. EXPORT FORMATTED RESULTS
# ================================================

def export_formatted_results(aggregated_df, sources):
    """
    Export publication-ready tables
    """
    # Create summary table
    summary_table = aggregated_df[['Gene']].copy()

    if 'Functional_tag' in aggregated_df.columns:
        summary_table['Functional_tag'] = aggregated_df['Functional_tag']

    # Add prevalence for each source
    for source in sources:
        summary_table[f'{source}_prev'] = aggregated_df[f'{source}_prevalence'].round(3)

    # Sort by overall prevalence
    summary_table['Overall_Prevalence'] = summary_table[[f'{s}_prev' for s in sources]].mean(axis=1)
    summary_table = summary_table.sort_values('Overall_Prevalence', ascending=False)

    # Save to Excel with formatting
    with pd.ExcelWriter('neuroactive_gene_results.xlsx', engine='openpyxl') as writer:
        # Main results
        summary_table.to_excel(writer, sheet_name='Gene_Prevalence', index=False)

        # Source statistics
        source_stats = []
        for source in sources:
            source_genes = summary_table[summary_table[f'{source}_prev'] > 0.5]
            source_stats.append({
                'Source': source,
                'Total_Genes': len(source_genes),
                'Mean_Prevalence': summary_table[f'{source}_prev'].mean(),
                'Top_Gene': source_genes.iloc[0]['Gene'] if len(source_genes) > 0 else 'None'
            })

        pd.DataFrame(source_stats).to_excel(writer, sheet_name='Source_Statistics', index=False)

        # Correlation matrix
        corr_data = []
        for source1 in sources:
            row = {'Source': source1}
            for source2 in sources:
                corr = np.corrcoef(summary_table[f'{source1}_prev'],
                                  summary_table[f'{source2}_prev'])[0, 1]
                row[source2] = round(corr, 3)
            corr_data.append(row)

        pd.DataFrame(corr_data).to_excel(writer, sheet_name='Correlation_Matrix', index=False)

    print("\n" + "="*60)
    print("EXPORT COMPLETE")
    print("="*60)
    print("\nGenerated files:")
    print("1. neuroactive_gene_heatmap.pdf - Main publication figure")
    print("2. neuroactive_gene_barcharts.pdf - Multi-panel analysis")
    print("3. interactive_heatmap.html - Interactive exploration")
    print("4. neuroactive_gene_results.xlsx - Complete results tables")
    print("5. statistical_summary.csv - Key statistics")
 

# Export all results
export_formatted_results(aggregated_df, sources)
