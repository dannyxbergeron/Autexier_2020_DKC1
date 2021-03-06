import pandas as pd

input_files = snakemake.input.bedgraphs
output_files = snakemake.output.clean_bg

colnames = ['chr', 'start', 'end', 'value']
chr = ['chr{}'.format(x+1) for x in range(22)] + ['chrY', 'chrX', 'chrM', 'chrMT']
for infile, outfile in zip(input_files, output_files):

    df = pd.read_csv(infile, sep='\t', names=colnames, skiprows=1, dtype={'chr': str})

    mask = df.chr.isin(chr)
    clean_df = df.loc[mask]

    clean_df.to_csv(outfile, sep='\t', index=False, header=False)


chrNameLength_df = pd.read_csv(snakemake.input.chrLength, sep='\t',
                               names=['chr', 'val'],
                               dtype={'chr': str})
chrNameLength_df['chr'] = 'chr' + chrNameLength_df['chr']
chrNameLength_df.to_csv(snakemake.output.new_chrLength,
                        sep='\t', index=False, header=False)
