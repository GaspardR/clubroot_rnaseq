#!/usr/bin/env python3

from pathlib import Path

import pandas as pd

ref_df = pd.read_csv(snakemake.input.map, "\t", names=["transcript", "gene"])

for file in snakemake.input.counts:
    name = Path(file).parent.name
    df = pd.read_csv(file, sep="\t")
    ref_df[name] = ref_df.transcript.map(dict(zip(df.target_id, df.est_counts)))

ref_df.drop(columns="transcript", inplace=True)
ref_df = ref_df.groupby("gene").sum()
ref_df.to_csv(snakemake.output.combined, header=True, index=False)
