import pandas as pd
import numpy as np

def tidy_data(df: pd.DataFrame) -> pd.DataFrame:
    to_drop = ['pValueAnnotation', 'beta', 'mappedGenes', 'ci', 'traitName', 'efoTraits', 'bgTraits', 'accessionId', 'locations', 'pubmedId', 'author', 'riskFrequency']
    tidy = df.copy()

    # For SNPs with an rsID, convert the names into standard chromosome format
    tidy['riskAllele'] = tidy['riskAllele'].mask(~tidy['riskAllele'].str.startswith("chr"), "chr" + tidy['locations'] + tidy['riskAllele'].str.slice(-2), axis=0)
    tidy = tidy[~tidy['riskAllele'].str.startswith("chr-")]

    # File from GWAS Catalog has '-' for missing values. Replace these with np.nan
    tidy.replace('-', np.nan, inplace=True)

    # Beta values are writen as 'unit/z-score increase/decrease'. Convert these to float
    tidy['beta'] = tidy['beta'].apply(lambda x: (-1 if "decrease" in x else 1) * float(x.split()[0]) if isinstance(x, str) else x)

    # Remove any extra text from orValue column
    tidy['orValue'] = tidy['orValue'].str.split(n=1, expand=True)[0]

    # Fill any missing values in orValue with the calculated value from beta
    tidy['orValue'] = tidy['orValue'].fillna(np.exp(tidy['beta']))

    # Drop any rows still missing orValue
    tidy.dropna(subset='orValue', inplace=True)

    # Convert column to float
    tidy['orValue'] = tidy['orValue'].astype(float)

    # Remove columns that are not needed
    tidy.drop(columns=to_drop, inplace=True)

    # Where multiple studies show the same association, only take the statistics for the study with the smallest p-value
    tidy = tidy.sort_values('pValue').drop_duplicates('riskAllele')
    return tidy

def main():
    pass

if __name__ == '__main__':
    main()