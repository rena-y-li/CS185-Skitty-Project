
# This file imports the scRNA-Seq data from the Onek1k cohort and creates pseudobulk matrices to be used during the CONTENT run

import scanpy as sc
import pandas as pd
import numpy as np
import os

# Read in data from directory

adata = sc.read_h5ad("expression/h5ad/onek1k.h5ad")

os.makedirs("real_pseudobulk", exist_ok=True)

# Extract all cell types

celltypes = adata.obs["predicted.celltype.l2"].unique()

for ct in celltypes:

    print(f"Processing: {ct}")

    sub = adata[adata.obs["predicted.celltype.l2"] == ct]

    # Extract unique donors
    
    donors = sub.obs["donor_id"].unique()

    pb_list = []

    for donor in donors:

        donor_sub = sub[sub.obs["donor_id"] == donor]

        # Average over cells in this donor
        
        expr = np.asarray(donor_sub.X.mean(axis=0)).ravel()

        pb_list.append(expr)

    # Rows = donors, Cols = genes
    
    pb = pd.DataFrame(
        pb_list,
        index=donors,
        columns=sub.var_names
    )

    safe_ct = ct.replace(" ", "_").replace("/", "_")

    pb.to_csv(
        f"real_pseudobulk/pseudobulk__{safe_ct}.tsv.gz",
        sep="\t",
        compression="gzip"
    )

