import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

GWAS_LINK = "https://www.ebi.ac.uk/gwas/api/v2/efotraits/{}/associations/download?includeBgTraits=false&includeChildTraits=true"


def tidy_associations(df: pd.DataFrame, keep=None) -> pd.DataFrame:
    if keep is None:
        keep = ['orValue', 'riskFrequency']

    to_drop = ['pValueAnnotation', 'beta', 'mappedGenes', 'ci', 'traitName', 'efoTraits', 'bgTraits', 'accessionId',
               'locations', 'pubmedId', 'author']
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
    to_drop = ['chromosome', 'base_pair_location', 'effect_allele', 'other_allele', 'effect_allele_frequency', 'beta', 'variant_id', 'ci_upper', 'ci_lower']
    tidy: pd.DataFrame = df[df['p_value'] < significance].copy()

    tidy['risk_allele'] = tidy['chromosome'].astype(str) + ":" + tidy['base_pair_location'].astype(str) + "-" + tidy['effect_allele'] + tidy['other_allele']

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
