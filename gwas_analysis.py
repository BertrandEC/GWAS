import pandas as pd
from gwas_analysis.stats import tidy_summary_stats

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