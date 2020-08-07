import pandas as pd

input_files = snakemake.input.tpm_files
output_counts = snakemake.output.counts
output_file = snakemake.output.merged

for i, file in enumerate(input_files):
    sample_name = file.replace('.tsv', '').split('/')[-1]

    if i == 0:
        master_df = pd.read_csv(file, sep='\t')
        counts_df = master_df.drop(columns=['cpm', 'tpm', 'gene_name'])
        counts_df.columns = ['gene_id', sample_name]

        tpm_df = master_df.drop(columns=['count', 'cpm'])
        tpm_df.columns = ['gene_id', 'gene_name', sample_name]
    else:
        df = pd.read_csv(file, sep='\t')
        counts_dict = dict(zip(df.gene_id, df.count))
        counts_dict[sample_name] = counts_dict.gene_id.map(counts_dict)

        tpm_dict = dict(zip(df.gene_id, df.tpm))
        tpm_dict[sample_name] = tpm_dict.gene_id.map(tpm_dict)

counts_df.to_csv(output_counts, sep='\t', index=False)
master_df.to_csv(output_file, sep='\t', index=False)
