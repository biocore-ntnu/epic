def get_number_of_islands_in_each_bin_per_sample(df):

    df = df[[c for c in df.columns if c.endswith("Island")]]
    return df
