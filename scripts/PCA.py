import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

input_file = snakemake.input.merged_tpm
output_tok = snakemake.output.tok

def load_data():

    df = pd.read_csv(input_file, sep='\t')
    df.index = df.gene_id
    df.drop(columns=['gene_id', 'gene_name'], axis=1, inplace=True)

    # # Filter for low abundance gene (does not seems to be good...)
    # print(len(df))
    # df['mean'] = df.mean(axis=1)
    # df = df.loc[df['mean'] > 50]
    # df.drop(columns=['mean'], inplace=True)
    # print(len(df))


    trans_df = df.transpose()

    return trans_df

def process_data(df):

    Y = df.index
    X = df.values
    print(X)
    # X = np.sqrt(X)
    X = StandardScaler().fit_transform(X)
    # print(Y)
    # print(X.shape)
    # print(X)

    return X, Y

def pca_analysis(X, Y):
    print(X)
    print(Y)

    conditions = {
        'mock' : '#e41a1c',
        'sidkc1' : '#377eb8',
        'sidkc1-WT' : '#4daf4a',
        'sidkc1-K446X' : '#984ea3',
        'sidkc1-K467R' : '#ff7f00',
    }

    pca = PCA(n_components=2)
    principalComponents= pca.fit_transform(X)
    print('Explained variation per principal component: {}'.format(pca.explained_variance_ratio_))
    # print(principalComponents)
    pca_df = pd.DataFrame(data=principalComponents,
                            columns = ['pc1',
                                        'pc2'], index=Y)
    pca_df.reset_index(inplace=True)
    pca_df[['condition', 'replicate']] = pca_df['index'].str.split('_', expand=True)
    print(pca_df)

    for s in conditions.keys():
        tmp = pca_df.loc[pca_df.condition == s]
        plt.scatter(tmp['pc1'], tmp['pc2'], label=s, s=100,
                    alpha=0.75, color=conditions[s], edgecolor='k')

    plt.xlabel('Principal Component 1 ({:.1f}%)'.format(pca.explained_variance_ratio_[0] * 100))
    plt.ylabel('Principal Component 2 ({:.1f}%)'.format(pca.explained_variance_ratio_[1] * 100))
    plt.title("Principal Component Analysis")
    plt.legend()
    plt.show()


def main():

    df = load_data()

    X, Y = process_data(df)

    pca_analysis(X, Y)



if __name__ == '__main__':
    main()
