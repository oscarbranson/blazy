from .io import run_phreeqc

def df_starchunking(df, chunksize, outputs, db):
    """Splits df into chunks, drops data of original df inplace"""
    while len(df):
        # Return df chunk
        yield df.iloc[:chunksize].copy(), outputs, db
        # Delete data in place because it is no longer needed
        df.drop(df.index[:chunksize], inplace=True)

