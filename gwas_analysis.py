import pandas as pd
from gwas_analysis.plots import plot_manhattan, plot_scatter_with_regression, plot_logistic_regression
from gwas_analysis.stats import tidy_summary_stats, logistic_regression_using_OR

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

    # plot_joint_manhattan(merged_df, 'risk_allele', 'p_value_diabetes', 'neg_log_p_value_arthritis')
    plot_manhattan(merged_df, 'p_value_diabetes', 'risk_allele')
    plot_scatter_with_regression(merged_df, col1='beta_diabetes', col2='beta_arthritis')

if __name__ == '__main__':
    main()