import pandas as pd
import numpy as np

def shorten_id(unirefid):
    # remove the
    return unirefid.split('_')[1]


def split_uid(mydf):
    """
    Some unirefids I will use uniprotkb to search. Other I will use uniref to search
    """

    mydf = mydf.iloc[1:, :] # remove UNMAPPED
    mydf['GeneShort'] = mydf['Gene Family'].map(shorten_id)

    uid_filter = ['UPI00' in uid for uid in mydf['GeneShort'].to_numpy()]

    mydf_2 = mydf.loc[uid_filter, :]
    mydf_1 = mydf.loc[np.logical_not(uid_filter), :]

    return mydf_1, mydf_2


if __name__ == "__main__":

    plaque_unirefs_df = pd.read_csv("unirefs/unique_uniref90_plaque.txt")
    plaque_unirefs_df_1, plaque_unirefs_df_2 = split_uid(plaque_unirefs_df)
    plaque_unirefs_df_1.to_csv("unirefs/uniref90_part1_plaque.csv", index=False)
    plaque_unirefs_df_2.to_csv("unirefs/uniref90_part2_plaque.csv", index=False)

    saliva_unirefs_df = pd.read_csv("unirefs/unique_uniref90_saliva.txt")
    saliva_unirefs_df_1, saliva_unirefs_df_2 = split_uid(saliva_unirefs_df)
    saliva_unirefs_df_1.to_csv("unirefs/uniref90_part1_saliva.csv", index=False)
    saliva_unirefs_df_2.to_csv("unirefs/uniref90_part2_saliva.csv", index=False)
