import pandas as pd
import matplotlib.pyplot as plt
import sys

"""
Script producing histograms from columns 'insert_size' and 'All_Reads.fr_count' contained in
Picard_insert_size_metrics.txt files. Each of these file must be in a separate folder named with the sample name. You
must change the sample_lst according to your sample names.
"""

input_files = snakemake.input.metrics
output_graph_files = snakemake.output.graphs
samples = snakemake.params.samples

def txt_to_df(file):
    """Extract useful data from txt file and return dataframe."""
    with open(file, 'r') as f:
        while not 'All_Reads.fr_count' in next(f):
            pass
        temp = [line.strip().split() for line in f]

        df = pd.DataFrame(temp, columns=['insert_size', 'all_read_counts'])
        df = df.dropna()
        for col in df.columns: df[col] = df[col].map(int)

        return df

def histogram(df, figure_output_path, sample):
    print(df)
    """Generate histogram from given dataframe."""
    ax = df.plot.bar(x='insert_size', y='all_read_counts', color='blue')
    ax.set_xlabel('Insert size', fontsize=12)
    ax.set_ylabel('Total number of inserts', fontsize=12)

    ax.set_xticklabels(df.insert_size, fontsize=10, rotation=0)
    ax.legend().remove()

    every_nth = 50
    for n, label in enumerate(ax.xaxis.get_ticklabels()):
        if n % every_nth != 0:
            label.set_visible(False)

    plt.title(sample)
    plt.tight_layout()
    # plt.show()
    plt.savefig(figure_output_path, dpi=300, bbox_inches='tight')


def main():
    """Iterate through sample list to generate a histogram for each sample (stored in the same folder than the txt
        file)"""
    for inp, out, s in zip(input_files, output_graph_files, samples):
        print(inp)
        print(out)
        print(s)
        df = txt_to_df(inp)
        histogram(df, out, s)


main()
