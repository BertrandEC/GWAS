import pandas as pd
from gwas_analysis.plots import plot_manhattan, plot_joint_manhattan, plot_scatter_with_regression, plot_logistic_regression, plot_prs_comparison
from gwas_analysis.stats import tidy_summary_stats, logistic_regression_using_OR, simulate_individual_prs_two_diseases, compute_prs_from_summary_stats

def main():
    pd.options.display.max_columns = None
    pd.options.display.max_rows = None

    try:
        merged_df = pd.read_csv('data/merged.tsv', sep='\t')
    except FileNotFoundError:
        t1dm = tidy_summary_stats(pd.read_csv('data/diabetes.tsv', sep='\t'), significance=1)
        ra = tidy_summary_stats(pd.read_csv('data/ra.tsv', sep='\t'), significance=1)

        merged_df = pd.merge(t1dm, ra, how='inner', on=['risk_allele', 'chromosome', 'base_pair_location', 'effect_allele', 'other_allele'], suffixes=("_diabetes", "_arthritis"))
        merged_df.to_csv('data/merged.tsv', index=False, sep='\t')

    plot_joint_manhattan(merged_df, 'neg_log_p_value_diabetes', 'neg_log_p_value_arthritis')
    #plot_manhattan(t1dm, ra, 'neg_log_p_value', 'neg_log_p_value')
    merged_df = merged_df[(merged_df['p_value_diabetes'] < 1e-5) | (merged_df['p_value_arthritis'] < 1e-5)]

    plot_scatter_with_regression(merged_df, col1='beta_diabetes', col2='beta_arthritis')
    plot_scatter_with_regression(merged_df, col1='neg_log_p_value_diabetes', col2='neg_log_p_value_arthritis')
    # plot_logistic_regression(logistic_regression_using_OR(merged_df, 'odds_ratio_diabetes', 'odds_ratio_arthritis'), merged_df['beta_diabetes'])


if __name__ == '__main__':
    main()