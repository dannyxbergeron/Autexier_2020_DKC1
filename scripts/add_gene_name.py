import os
import pandas as pd

from snakemake.shell import shell

input_dir = snakemake.params.genes_dir
inputs = [
    x
    for x in os.listdir(os.path.join(input_dir))
    if x.endswith('csv')
    and 'analysed_' not in x
]
corresponding = snakemake.params.corresponding

tpms = snakemake.params.tpms
tpm_df = pd.read_csv(tpms, sep='\t')
tpm_df.set_index('gene',inplace=True)
tpm_df.columns = [x[:-2] for x in tpm_df]
tpm_df = tpm_df.transpose()
tpm_df.reset_index(inplace=True)
tpm_df = tpm_df.groupby('index').mean().transpose()

corr_df = pd.read_csv(corresponding, sep='\t', names=['gene_id', 'gene_name'])
corr_dict = dict(zip(corr_df.gene_id, corr_df.gene_name))

for in_file in inputs:
    file = os.path.join(input_dir, in_file)
    df = pd.read_csv(file)
    df.rename(columns={'Unnamed: 0': 'gene_id'}, inplace=True)
    df.drop(columns=['baseMean', 'lfcSE', 'stat', 'pvalue'], inplace=True)
    df['gene_name'] = df.gene_id.map(corr_dict)

    df.dropna(axis=0, inplace=True)
    df.sort_values('padj', ascending=True, inplace=True)
    df = df[[df.columns[0], df.columns[-1], *df.columns[1:-1]]]

    df['Lipo'] = df.gene_id.map(dict(zip(tpm_df.index, tpm_df.Lipo)))
    df['KD'] = df.gene_id.map(dict(zip(tpm_df.index, tpm_df.NOP58)))

    df = df.loc[(df.padj < 0.05) & ((df.log2FoldChange > 1) | (df.log2FoldChange < -1))]

    df.to_csv(os.path.join(input_dir, 'filtered_analysed_' + in_file), sep='\t', index=False)

shell('touch {}'.format(snakemake.output.tok))
