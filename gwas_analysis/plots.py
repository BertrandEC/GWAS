import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from stats import manual_pearson_correlation

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


def plot_manhattan(df, pval_col, snp_col):
    # Sort by SNP
    df_sorted = df.sort_values(by=snp_col)
    # Plot the Manhattan plot
    plt.figure(figsize=(10, 6))
    plt.scatter(df_sorted[snp_col], -np.log10(df_sorted[pval_col]), c='blue', s=10)
    plt.xlabel('SNP')
    plt.ylabel('-log10(p-value)')
    plt.title('Manhattan Plot for SNP Significance')
    plt.show()


def plot_joint_manhattan(df, snp_pos_col, neg_log_pval_col1, neg_log_pval_col2):
    """Plots a joint Manhattan plot comparing SNPs possitions across two diseases using -log10(p-values). """

    plt.figure(figsize=(12, 6))
    plt.scatter(df[snp_pos_col], df[neg_log_pval_col1], color='blue', s=10, label='Disease 1', alpha=0.6)
    plt.scatter(df[snp_pos_col], df[neg_log_pval_col2], color='red', s=10, label='Disease 2', alpha=0.6)

    # Add labels
    plt.xlabel('SNP positions')
    plt.ylabel('-log10(p-value)')
    plt.title('Joint Manhattan Plot Comparing Disease 1 and Disease 2')
    plt.legend()
    plt.show()


def plot_scatter_with_regression(df, col1, col2):
    """For visualizing relationships between two variables"""
    correlation = manual_pearson_correlation(df, col1, col2)

    plt.figure(figsize=(8, 6))
    sns.regplot(x=col1, y=col2, data=df, line_kws={"color": "red"}, scatter_kws={"s": 50})
    plt.xlabel(col1)
    plt.ylabel(col2)
    plt.title(f'Scatter plot with Regression line: {col1} vs {col2}\n Pearson r = {correlation:.3f}')
    plt.show()


def plot_logistic_regression(model, x_vals):
    """visualize the relationship between log(OR1) and the probability of SNP significance in Disease 2.
    A common way to do this is by plotting the sigmoid function."""
    # Generate the log-odds from model predictions (β₀ + β₁ * log(OR1))
    log_odds = model.params[0] + model.params[1] * np.log(x_vals)
    # Convert log-odds to probability using the sigmoid function
    probabilities = 1 / (1 + np.exp(-log_odds))

    plt.figure(figsize=(8, 6))
    plt.plot(x_vals, probabilities, label='Logistic Regression Curve', color='blue')
    plt.xlabel('log(OR1)')
    plt.ylabel('Probability of SNP significance in Disease 2')
    plt.title('Logistic Regression: Probability of SNP significance')
    plt.grid(True)
    plt.show()