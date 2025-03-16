import numpy as np
import pandas as pd
import statsmodels.api as sm
import scipy.stats as stats


def tidy_summary_stats(df: pd.DataFrame, significance=1e-2, rsid=False) -> pd.DataFrame:
    to_keep = ['risk_allele', 'chromosome', 'base_pair_location', 'p_value', 'neg_log_p_value', 'beta', 'odds_ratio', 'z_score', 'effect_allele_frequency', 'effect_allele', 'other_allele']
    tidy: pd.DataFrame = df[df['p_value'] < significance].copy()

    rsid_col = next(filter(lambda c: c in tidy, ['variant_id', 'rsid']), None)
    if rsid or rsid_col is not None:
        tidy['risk_allele'] = tidy[rsid_col]
    else:
        tidy['risk_allele'] = tidy['chromosome'].astype(str) + ":" + tidy['base_pair_location'].astype(str) + ":" + tidy['effect_allele'] + ':' + tidy['other_allele']

    tidy['neg_log_p_value'] = -np.log10(tidy['p_value'])
    if 'odds_ratio' in tidy:
        if not 'beta' in tidy:
            tidy['odds_ratio'] = np.nan
        tidy['beta'] = tidy['beta'].fillna(np.log(tidy['odds_ratio']))
    if 'beta' in tidy:
        if not 'odds_ratio' in tidy:
            tidy['odds_ratio'] = np.nan
        tidy['odds_ratio'] = tidy['odds_ratio'].fillna(np.exp(tidy['beta']))

    if 'standard_error' in tidy:
        tidy['z_score'] = tidy['beta'] / tidy['standard_error']
    else:
        tidy['z_score'] = np.nan
    tidy['z_score_p'] = np.sqrt(stats.chi2.isf(tidy['p_value'], 1))
    tidy['z_score'] = tidy['z_score'].fillna(tidy['z_score_p'])

    # Remove columns that are not needed
    tidy = tidy[to_keep]
    return tidy
 
def manual_pearson_correlation(df, x, y):
    """Can use to check linear relationships between odds ratios or p-values.
    correlation tells you "how much" two variables are related."""
    if isinstance(x, str):
        x = df[x]
        y = df[y]

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
    return manual_pearson_correlation(df, x, y)


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