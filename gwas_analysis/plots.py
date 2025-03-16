import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from .stats import manual_pearson_correlation, manual_spearman_correlation, manual_linear_regression

def plot_prs_comparison(prs1, prs2):
    """Plots the PRS distributions for two diseases to compare genetic risk."""
    plt.figure(figsize=(8, 6))
    sns.histplot(prs1, kde=True, color='blue', label="Disease 1", bins=50)
    sns.histplot(prs2, kde=True, color='red', label="Disease 2", bins=50)
    plt.xlabel('Polygenic Risk Score (PRS)')
    plt.ylabel('Density')
    plt.title('PRS Distribution Comparison')
    plt.legend()
    plt.show()


def plot_manhattan(df: pd.DataFrame, df2: pd.DataFrame, neg_log_pval_col1: str, neg_log_pval_col2: str):
    df_sorted  = df.sort_values(['chromosome','base_pair_location'])
    df2_sorted  = df2.sort_values(['chromosome','base_pair_location'])
    df_sorted.reset_index(inplace=True, drop=True)
    df2_sorted.reset_index(inplace=True, drop=True)
    df_sorted['i'] = df_sorted.index
    df2_sorted['i'] = df2_sorted.index + len(df_sorted)

    combined_df = pd.concat([df_sorted, df2_sorted])

    combined_df['dataset'] = ['Dataset 1'] * len(df_sorted) + ['Dataset 2'] * len(df2_sorted)

    # Create the Manhattan plot using sns.relplot
    plot = sns.relplot(
        data=combined_df,
        x='i',
        y=neg_log_pval_col1,  # Use the first dataset's p-value column
        hue='chromosome',
        style='dataset',  # Distinguish datasets by marker style
        palette='bright',
        aspect=3.7,
        height=6,
        legend=False
    )
    chrom_df = combined_df.groupby('chromosome')['i'].median()
    plot.ax.set_xticks(chrom_df)
    plot.ax.set_xticklabels(chrom_df.index)
    plot.ax.set_xlabel('Chromosome')
    plot.ax.set_ylabel('-log10(P-value)')
    plot.fig.suptitle('Joint Manhattan Plot')

    # Add genome-wide significance threshold
    plot.ax.axhline(-np.log10(5e-8), color='red', linestyle='--', label='Genome-wide significance')
    plot.ax.legend()

    plt.show()



def plot_joint_manhattan(df: pd.DataFrame, neg_log_pval_col1: str, neg_log_pval_col2: str):
    """Plots a joint Manhattan plot comparing SNPs possitions across two diseases using -log10(p-values). """

    df = df.sort_values(['chromosome', 'base_pair_location'])

    # Compute cumulative base pair positions for correct spacing across chromosomes
    chrom_offsets = {}
    offset = 0
    cumulative_positions = []

    for chrom in sorted(df['chromosome'].unique()):
        chrom_data = df[df['chromosome'] == chrom]
        chrom_offsets[chrom] = offset
        cumulative_positions.extend(chrom_data['base_pair_location'] + offset)
        offset += chrom_data['base_pair_location'].max()

    df['cumulative_pos'] = cumulative_positions

    # Plot the data
    plt.figure(figsize=(10, 5))
    plt.scatter(df['cumulative_pos'], df[neg_log_pval_col1], color='blue', s=10, label='Disease 1', alpha=0.6)
    plt.scatter(df['cumulative_pos'], df[neg_log_pval_col2], color='red', s=10, label='Disease 2', alpha=0.6)

    # Add chromosome labels at the middle of each chromosome section
    chrom_labels = {}
    for chrom in chrom_offsets:
        chrom_start = chrom_offsets[chrom]
        chrom_end = df[df['chromosome'] == chrom]['cumulative_pos'].max()
        chrom_labels[chrom] = (chrom_start + chrom_end) / 2
        plt.axvline(x=chrom_start, color='black', linestyle='--', linewidth=1)  # Vertical divider

    # Set x-axis ticks at chromosome midpoints
    plt.xticks(list(chrom_labels.values()), list(chrom_labels.keys()))
    plt.xlabel('Chromosome')
    plt.ylabel('-log10(p-value)')
    plt.title('Joint Manhattan Plot Comparing Disease 1 and Disease 2')

    plt.legend()
    plt.tight_layout()
    plt.show()



def plot_scatter_with_regression(df: pd.DataFrame, col1: str, col2: str):
    """For visualizing relationships between two variables"""
    correlation1 = manual_pearson_correlation(df, col1, col2)
    correlation2 = manual_spearman_correlation(df, col1, col2)

    slope, intercept = manual_linear_regression(df, col1, col2)
    

    plt.figure(figsize=(8, 6))
    sns.regplot(x=col1, y=col2, data=df, line_kws={"color": "red"}, scatter_kws={"s": 50})

    equation_text = f'y = {slope:.3f}x + {intercept:.3f}'
    plt.text(x=df[col1].min(), y=df[col2].max(), s=equation_text, fontsize=12, color="black", bbox=dict(facecolor='white', alpha=0.5))

    plt.xlabel(col1)
    plt.ylabel(col2)
    plt.title(f'Scatter plot with Regression line: {col1} vs {col2}\n Pearson r = {correlation1:.3f}, Spearman s = {correlation2:.3f}')
    plt.show()


def plot_logistic_regression(model, x_vals):
    """visualize the relationship between log(OR1) and the probability of SNP significance in Disease 2.
    A common way to do this is by plotting the sigmoid function."""
    # Generate the log-odds from model predictions (β₀ + β₁ * log(OR1))
    log_odds = model.params[0] + model.params[1] * x_vals
    # Convert log-odds to probability using the sigmoid function
    probabilities = 1 / (1 + np.exp(-log_odds))

    plt.figure(figsize=(8, 6))
    plt.plot(x_vals, probabilities, label='Logistic Regression Curve', color='blue')
    plt.xlabel('log(OR1)')
    plt.ylabel('Probability of SNP significance in Disease 2')
    plt.title('Logistic Regression: Probability of SNP significance')
    plt.grid(True)
    plt.show()