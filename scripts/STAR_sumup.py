import pandas as pd

INPUTS = snakemake.input.log_files
OUTPUT = snakemake.output.merged

def read_log_final_out(file,name):
    log_final_out_df=pd.read_csv(file,sep='|',header=None, names=['info',name])
    log_final_out_df['info']=log_final_out_df['info'].str.strip()
    log_final_out_df[name]=log_final_out_df[name].str.strip()
    return log_final_out_df


def main():

    for i,file in enumerate(INPUTS):
        name=file.split('/')[-2]
        log_final_out_df=read_log_final_out(file,name)
        if i==0:
            log_final_out_df_final=log_final_out_df.copy(deep=True)
        else:
            log_final_out_df_final=pd.merge(log_final_out_df_final,
                                                log_final_out_df,
                                                how='left',
                                                on='info')
    log_final_out_df_final.to_csv(OUTPUT,
                  index=False, header=True, sep=',')




main()
