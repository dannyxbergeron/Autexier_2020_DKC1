import numpy as np
import pandas as pd

input_files = snakemake.input.samtools_idx
output_file = snakemake.output.merged

colnames = ['length', 'nb_reads', 'unknown']
for i, file in enumerate(input_files):
    name = '_'.join(file.split('/')[-1].split('_')[:2])
    print(i, name)

    df = pd.read_csv(file, sep='\t', names=colnames, dtype={'chr': str})

    if i == 0:
        master_df = df[['length', 'nb_reads']]
        master_df.columns = ['length', name]
    else:
        reads_dict = dict(zip(df.index, df.nb_reads))
        master_df[name] = master_df.index.map(reads_dict)

master_df.reset_index(inplace=True)
master_df.rename({'index': 'chr'}, inplace=True, axis=1)

# Adding the total reads count at the end of the dataframe
to_append = ['Total_reads', np.nan] + [master_df[col].sum() for col in master_df.columns[2:]]
total_reads = pd.Series(to_append, index=master_df.columns)
master_df = master_df.append(total_reads, ignore_index=True)

# Adding the sum of ERCC reads for each sample
ERCC = master_df.loc[master_df.chr.str.startswith('ERCC')]
ERCC_sum = ['ERCC_reads', np.nan] + [ERCC[col].sum() for col in ERCC.columns[2:]]
ERCC_sum = pd.Series(ERCC_sum, index=master_df.columns)
master_df = master_df.append(ERCC_sum, ignore_index=True)

# Adding the pourcent of ERCC reads
wanted_cols = list(master_df.columns[2:])
t_reads, ERCC_reads = master_df[wanted_cols].loc[master_df.chr.isin(['Total_reads', 'ERCC_reads'])].values
ERCC_percent = ERCC_reads / t_reads * 100
ERCC_percent = pd.Series(['ERCC_percent', np.nan] + list(ERCC_percent), index=master_df.columns)
master_df = master_df.append(ERCC_percent, ignore_index=True)

master_df.to_csv(output_file, sep='\t', index=False)
