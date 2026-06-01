# This script transposes the pseudobulk matrices to the correct format for CONTENT (wrong orientation in make_pseudobulk.py)

import pandas as pd
import glob
import os

os.makedirs("fixed_pseudobulk", exist_ok=True)

for f in glob.glob("real_pseudobulk/*.tsv.gz"):

    df = pd.read_csv(f, sep="\t", index_col=0)

    # Transpose so that the genes become rows
    df = df.T

    # The gene IDs become the first column
    df.insert(0, "gene_id", df.index)

    # Reset row numbering
    df = df.reset_index(drop=True)

    out = os.path.basename(f)

    # Export fixed files
    df.to_csv(
        f"fixed_pseudobulk/{out}",
        sep="\t",
        index=False,
        compression="gzip"
    )

    print(out, df.shape)
