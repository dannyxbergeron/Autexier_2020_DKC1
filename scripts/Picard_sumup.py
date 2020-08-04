import pandas as pd

INPUTS = snakemake.input.metrics
OUTPUT = snakemake.output.merged

pd.set_option('max_rows', 100)
pd.set_option('max_columns', 80)
pd.set_option('display.width', 400)

def read_Picard_file(file,name):
    with open(file,'r') as p:
        Picard_content=p.readlines()
    start_index=[i for i in range(len(Picard_content)) if Picard_content[i].startswith('MEDIAN_INSERT_SIZE')==True][0]
    columns=list(Picard_content[start_index].rstrip().split())[:6]
    values=list(Picard_content[start_index+1].rstrip().split())[:6]
    Picard_df=pd.DataFrame(data={'info':columns,name:values})
    return Picard_df


def main():

    for i,file in enumerate(INPUTS):
        name=file.split('/')[-2]
        Picard_df=read_Picard_file(file,name)
        if i==0:
            Picard_df_final=Picard_df.copy(deep=True)
        else:
            Picard_df_final=pd.merge(Picard_df_final,Picard_df,how='left',on='info')

    cols = list(Picard_df_final)
    cols.insert(0, cols.pop(cols.index('info')))
    Picard_df_final = Picard_df_final.ix[:, cols]
    print(Picard_df_final[:5])

    Picard_df_final.to_csv(OUTPUT,
                  index=False, header=True, sep=',')


main()
