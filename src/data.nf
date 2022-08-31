nextflow.enable.dsl=2

// Convert Seurat object to h5ad file
process prepare_rna_data {
    publishDir 'result'

    output:
    file("RNA.h5ad")

    """
    #!/usr/bin/env Rscript
    library("anndata")
    library("Seurat")
    data1 <- readRDS("/home/kaizhang/data/project/mouse_brain/data/RNA/excitatory_IDonly_mergedAges.RDS")
    adata1 <- AnnData(
        X = GetAssayData(data1),
        var = data.frame(my.cell.type = data1\$my.cell.type), 
    )
    data2 <- readRDS("/home/kaizhang/data/project/mouse_brain/data/RNA/inhibitory_IDonly_mergedAges.RDS") 
    adata2 <- AnnData(
        X = GetAssayData(data2),
        var = data.frame(my.cell.type = data2\$my.cell.type), 
    )
    adata <- concat(c(adata1\$T, adata2\$T))
    adata\$write_h5ad("RNA.h5ad")
    """
}

// Convert metadata
process prepare_meta_data {
    publishDir 'result/data'

    output:
    file("metadata.tsv.gz")

    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri
    readRDS = robjects.r['readRDS']
    meta1 = readRDS("${metadata_excitatory}")
    meta1 = pandas2ri.rpy2py_dataframe(meta1)
    meta2 = readRDS("${metadata_inhibitory}")
    meta2 = pandas2ri.rpy2py_dataframe(meta2)
    df = pd.concat([meta1, meta2])
    df['orig.ident'] = df.index
    df.index = df['sample'] + "+" + np.array(list(map(lambda x: x.split("_")[-1], df.index)))
    df.to_csv("metadata.tsv.gz", index_label="ID", sep='\t')
    """
}