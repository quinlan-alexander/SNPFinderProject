#!/usr/bin/env python3
# Name: Quinlan Alexander (qalexand)
# Group Members: None

import pandas as pd

# /Users/ba25714/PycharmProjects/final_project

samples_df = pd.read_csv('/Users/ba25714/PycharmProjects/final_project/sample_results.csv')

score_counts_by_gene_df = pd.DataFrame(columns=['gene', 'warning', 'caution', 'ok'])

for gene in samples_df['gene'].unique():
    warn = samples_df.loc[(samples_df['gene'] == gene) & samples_df['score'] <= .05].agg(['count'])['score'][0]
    caution = samples_df.loc[(samples_df['gene'] == gene) & (samples_df['score'] > .05) & (samples_df['score'] <= 0.5)].agg(['count'])['score'][0]
    ok = samples_df.loc[(samples_df['gene'] == gene) & (samples_df['score'] > .50) & (samples_df['score'] <= 1.0)].agg(['count'])['score'][0]

    score_counts_by_gene_df = score_counts_by_gene_df.append({'gene': gene, 'warning': warn, 'caution': caution, 'ok': ok}, ignore_index=True)

score_counts_by_gene_df.style

# another way using indexes

def sum_warn_by(gene):
    return samples_df.loc[(samples_df['gene'] == gene) & samples_df['score'] <= .05].agg(['count'])['score'][0]

def sum_caution_by(gene):
    return samples_df.loc[(samples_df['gene'] == gene) & (samples_df['score'] > .05) & (samples_df['score'] <= 0.5)].agg(['count'])['score'][0]

def sum_ok_by(gene):
    return samples_df.loc[(samples_df['gene'] == gene) & (samples_df['score'] > .50) & (samples_df['score'] <= 1.0)].agg(
        ['count'])['score'][0]

def sum_sift_categories_by():
    bar = [pd.Series([gene, sum_warn_by(gene), sum_caution_by(gene), sum_ok_by(gene)]).reshape(1, 4) for gene in
                            samples_df['gene'].unique()
    foo = pd.DataFrame([pd.Series([gene, sum_warn_by(gene), sum_caution_by(gene), sum_ok_by(gene)]).reshape(1, 4) for gene in
                            samples_df['gene'].unique()])
    return foo

#sum_sift_categories_by()


score_counts_by_gene_df = pd.DataFrame([pd.Series([gene, sum_warn_by(gene), sum_caution_by(gene), sum_ok_by(gene)]).reshape(1, 4) for gene in
                            samples_df['gene'].unique()])

for gene in samples_df['gene'].unique():
    score_counts_by_gene_df = score_counts_by_gene_df.append({'gene': gene, 'warning': sum_warn_by(gene), 'caution': sum_caution_by(gene), 'ok': sum_ok_by(gene)}, ignore_index=True)


score_counts_by_gene_df.style