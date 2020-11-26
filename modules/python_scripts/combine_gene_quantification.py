#!/usr/bin/env python3

from pathlib import Path

import pandas as pd

files_df = [
    (
        pd.read_csv(
            file,
            sep='\t',
            header=0,
            index_col=False
        ),
        Path(file).name[:-4]
    )
    for file in snakemake.input.counts
]

ref_df = files_df[0][0]
data_df = ref_df[['gene_id']]
for i, file in enumerate(files_df):
    df, name = file
    df.columns = map(str.lower, df.columns)
    data_df[name] = data_df.gene_id.map(
        dict(zip(df.gene_id, df['count']))
    )
data_df = data_df.rename(columns={'gene_id':'gene'})
data_df.to_csv(snakemake.output.combined, sep=',', header=True, index=False)
