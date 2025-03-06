import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import seaborn as sns

GWAS_LINK = "https://www.ebi.ac.uk/gwas/api/v2/efotraits/{}/associations/download?includeBgTraits=false&includeChildTraits=true"


def tidy_associations(df: pd.DataFrame, keep=None) -> pd.DataFrame:
    if keep is None:
        keep = ['orValue', 'riskFrequency']

    to_drop = ['pValueAnnotation', 'beta', 'mappedGenes', 'ci', 'traitName', 'efoTraits', 'bgTraits', 'accessionId',
               'locations', 'pubmedId', 'author', 'pValue']
    tidy = df.copy()

    # For SNPs with an rsID, convert the names into standard chromosome format
    tidy['riskAllele'] = tidy['riskAllele'].mask(~tidy['riskAllele'].str.startswith("chr"),
                                                 "chr" + tidy['locations'] + tidy['riskAllele'].str.slice(-2), axis=0)
    tidy = tidy[~tidy['riskAllele'].str.startswith("chr-")]

    # File from GWAS Catalog has '-' for missing values. Replace these with np.nan
    tidy.replace('-', np.nan, inplace=True)
    tidy.replace('NR', np.nan, inplace=True)

    # Beta values are writen as 'x unit/z-score increase/decrease'. Convert these to float
    tidy['beta'] = tidy['beta'].apply(
        lambda x: (-1 if "decrease" in x else 1) * float(x.split()[0]) if isinstance(x, str) else x)

    # Remove any extra text from orValue column
    tidy['orValue'] = tidy['orValue'].str.split(n=1, expand=True)[0]

    # Fill any missing values in orValue with the calculated value from beta (or = exp(beta))
    tidy['orValue'] = tidy['orValue'].fillna(np.exp(tidy['beta']))

    # Drop any rows still missing orValue
    tidy.dropna(subset=keep, how='all', inplace=True)

    # Convert orValue column to float
    tidy['orValue'] = tidy['orValue'].astype(float)

    tidy['neg_log_p_value'] = tidy['pValue'].apply(lambda x: -np.log10(x))

    # Remove columns that are not needed
    tidy.drop(columns=to_drop, inplace=True)

    # Where multiple studies show the same association, only take the statistics for the study with the smallest p-value
    tidy = tidy.sort_values('pValue').drop_duplicates('riskAllele')

    tidy.rename(columns={
        'orValue': 'odds_ratio',
        'riskAllele': 'risk_allele',
    }, inplace=True)

    tidy.set_index('risk_allele', inplace=True)
    return tidy


def tidy_summary_stats(df: pd.DataFrame, significance=1e-2) -> pd.DataFrame:
    to_drop = ['p_value', 'chromosome', 'base_pair_location', 'effect_allele', 'other_allele', 'effect_allele_frequency', 'beta', 'variant_id', 'ci_upper', 'ci_lower']
    tidy: pd.DataFrame = df[df['p_value'] < significance].copy()

    tidy['risk_allele'] = tidy['chromosome'].astype(str) + ":" + tidy['base_pair_location'].astype(str) + "-" + tidy['effect_allele'] + tidy['other_allele']

    tidy['neg_log_p_value'] = tidy['p_value'].apply(lambda x: -np.log10(x))

    # Remove columns that are not needed
    tidy.drop(columns=to_drop, inplace=True, errors='ignore')
    tidy.set_index('risk_allele', inplace=True)

    return tidy


def dataframe_from_gwas(efo_id: str) -> pd.DataFrame:
    return tidy_associations(pd.read_csv(GWAS_LINK.format(efo_id), sep='\t'))


def main():
    pd.options.display.max_columns = None
    pd.options.display.max_rows = None

    t1dm = tidy_summary_stats(pd.read_csv('data/diabetes.tsv', sep='\t'))
    ra = tidy_summary_stats(pd.read_csv('data/ra.tsv', sep='\t').loc[:, 'variant_id':'beta'])

    merged_df = pd.merge(t1dm, ra, how='inner', on='risk_allele', suffixes=("_diabetes", "_arthritis"))
    merged_df.to_csv('data/merged.tsv', sep='\t')

    print(merged_df)


if __name__ == '__main__':
    main()


def manual_pearson_correlation(df, col1, col2):
    """Can use to check linear relationships between odds ratios or p-values.
    correlation tells you "how much" two variables are related."""
    x = df[col1]
    y = df[col2]

    x_mean = np.mean(x)
    y_mean = np.mean(y)

    numerator = np.sum((x - x_mean) * (y - y_mean))
    denominator = np.sqrt(np.sum((x - x_mean) ** 2) * np.sum((y - y_mean) ** 2))

    return numerator / denominator if denominator != 0 else np.nan
""" orValue_disease1 vs orValue_disease2: Measures shared genetic effects (should use log)
    neg_log_p_disease1 vs neg_log_p_disease2: Determines SNP significance overlap
    riskFrequency_disease1 vs riskFrequency_disease2: Examines shared allele frequencies
    orValue_disease1 vs riskFrequency_disease2: Tests if high-impact SNPs in one disease are more common in another can do vice vesra as well"""


def manual_spearman_correlation(df, col1, col2):
    """Computes Spearman correlation which converts data into ranks before computing Pearson correlation."""
    # Convert values to ranks
    x = df[col1].rank()
    y = df[col2].rank()

    # Compute Pearson correlation on ranked values
    return manual_pearson_correlation(x, y)


def manual_linear_regression(df, col1, col2):
    """Calculates the slope and intercept of the best-fit line for two columns in a dataframe.
    Linear Regression tells you if high OR SNPs in Disease 1 have high OR in Disease 2"""
    x = df[col2]
    y = df[col1]

    x_mean = np.mean(x)
    y_mean = np.mean(y)

    # Compute slope
    numerator = np.sum((x - x_mean) * (y - y_mean))
    denominator = np.sum((x - x_mean) ** 2)
    slope = numerator / denominator

    # Compute intercept
    intercept = y_mean - slope * x_mean

    return slope, intercept
""" orValue_disease1 vs orValue_disease2:
    If the slope is positive and close to 1, there is a strong shared genetic effect between the two diseases.(should use log)
    neg_log_p_disease1 vs neg_log_p_disease2:
    If there is a strong linear relationship, it suggests that significant SNPs in one disease tend to be significant in the other.
    riskFrequency_disease1 vs riskFrequency_disease2:
    If the slope is positive and close to 1, it suggests that the allele frequencies of high-impact SNPs are similar between the two diseases.
    orValue_disease1 vs riskFrequency_disease2:
    If the slope is positive and close to 1, it suggests that high-impact SNPs in one disease are more common in another. Can do vice versa as well."""


def logistic_regression_using_OR(df, or_col_1, or_col_2):
    """
    Estimate a logistic regression model between diseases using log(OR_2) = β_0 + β_1 * log(OR_1).
    Returning model (sm.Logit): Fitted logistic regression model.
    Logistic regression tells you if high OR SNPs in Disease 1 predict SNP significance in Disease 2.
    """

    X = np.log(df[or_col_1])  # Predictor: log(OR) of disease 1
    y = (df[or_col_2] > 1).astype(int)  # Response: Binary outcome (1 if associated, 0 if not)
    # could use a defined P-value threshold

    X = sm.add_constant(X)  # Add intercept
    model = sm.Logit(y, X).fit()

    return model
"""The coefficient (beta1) tells us how much log(OR) of disease 1 predicts disease 2 association.
If beta1 > 0, a higher OR in disease 1 increases the probability of the SNP being significant in disease 2.
If beta1 < 0, a higher OR in disease 1 reduces the probability of association in disease 2.
You can interpret the results for individual SNPs using the logistic regression coefficients:

Log-odds for each SNP: Multiply the SNP's log(OR1) value by the regression coefficient for log(OR1).
Probability: You can compute the probability of the SNP being significant in Disease 2 given its log(OR1):
P = 1 / (1 + exp(-β0 - β1 * log(OR1)))"""


def compute_prs_from_summary_stats(df, or_col, freq_col=None):
    """
    Computes Polygenic Risk Scores (PRS) using only GWAS summary statistics (effect size and allele frequencies).
    Returns the aggregate PRS for the given dataset.
    It calculates the Polygenic Risk Score for a given population or dataset based on available GWAS summary statistics.
    """
    # Calculate the PRS as the weighted sum of SNP effect sizes and allele frequencies
    freq = df[freq_col] if freq_col is not None else 0.5
    Beta = np.log(df[or_col])
    prs = np.sum(Beta * freq)
    prs_normalized = prs / len(df)

    return prs_normalized
"""We could also instead use seperate individual datasets of genotypes to calculate PRS individually.
If we are not given the frequencies we can assume they are 0.5 for each SNP."""


"""(categorical association) could also make a function that looks at the pre p_value filtered merged data frame 
and calculates the probability that a SNP is significant in Disease 2, given that it is already significant in Disease 1."""

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


def simulate_individual_prs_two_diseases(df, or_col1, freq_col1, or_col2, freq_col2, num_individuals=1000):
    """Simulates individual-level PRS for two diseases using GWAS summary statistics."""
    # Compute the effect sizes as the log of the odds ratios
    Beta1 = np.log(df[or_col1])  # Effect sizes for Disease 1
    Beta2 = np.log(df[or_col2])  # Effect sizes for Disease 2

    # Default allele frequency is 0.5 if not provided
    freq1 = df[freq_col1] if freq_col1 is not None else 0.5
    freq2 = df[freq_col2] if freq_col2 is not None else 0.5

    # Simulate genotypes using the binomial distribution based on Hardy-Weinberg equilibrium.
    # Simulate genotypes (0, 1, or 2 copies of risk allele) for num_individuals
    # Using binomial distribution: number of successes (alleles) in 2 trials, with probability `freq`
    genotypes1 = np.random.binomial(2, freq1.values[:, np.newaxis], (len(df), num_individuals))
    genotypes2 = np.random.binomial(2, freq2.values[:, np.newaxis], (len(df), num_individuals))

    # Compute the Polygenic Risk Scores (PRS) for each individual
    prs_values1 = np.dot(genotypes1.T, Beta1)  # PRS for Disease 1
    prs_values2 = np.dot(genotypes2.T, Beta2)  # PRS for Disease 2

    return prs_values1, prs_values2


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
