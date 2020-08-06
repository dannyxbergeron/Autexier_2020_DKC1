import numpy as np
import pandas as pd

tpm_file = snakemake.input.merged_tpm
ribo_file = snakemake.input.ribosomal_genes
output_file = snakemake.output.rRNA_values

df = pd.read_csv(tpm_file, sep='\t')
ribo_ids_df = pd.read_csv(ribo_file, sep='\t')
ribo_ids = list(ribo_ids_df.gene_id)

all_df = df.drop(columns=['gene_id', 'gene_name']).copy(deep=True)
all_sum = all_df.sum(axis=0)

ribo_df = df.loc[df.gene_id.isin(ribo_ids)]
ribo_df = ribo_df.drop(columns=['gene_id', 'gene_name']).copy(deep=True)
ribo_sum = ribo_df.sum(axis=0)
print(ribo_sum)

new_df = pd.DataFrame(all_sum, columns=['Total_TPM'])
new_df['rRNA_TPM'] = ribo_sum
new_df['%rRNA'] = new_df['rRNA_TPM'] / new_df['Total_TPM'] * 100

new_df.to_csv(output_file, sep='\t')
