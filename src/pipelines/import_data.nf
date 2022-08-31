/* Import and preprocess raw data files.
*/

nextflow.enable.dsl=2

params.baseOutputDir = 'result/imported_data/'

////////////////////////////////////////////////////////////////////////////////
// RNA data import and preprocessing
////////////////////////////////////////////////////////////////////////////////

// Convert Seurat object to h5ad file
process import_rna {
    publishDir "${params.baseOutputDir}/RNA/"

    // inputs are rds files provided by Harris
    input:
    file(rds_file_excitatory)
    file(rds_file_inhibitory)

    output:
    file("RNA.h5ad")

    """
    #!/usr/bin/env Rscript
    library("anndata")
    library("Seurat")
    data_excitatory <- readRDS("${rds_file_excitatory}")
    data_excitatory <- AnnData(
        X = GetAssayData(data_excitatory),
        var = data.frame(my.cell.type = data_excitatory\$my.cell.type), 
    )
    data_inhibitory <- readRDS("${rds_file_inhibitory}")
    data_inhibitory <- AnnData(
        X = GetAssayData(data_inhibitory),
        var = data.frame(my.cell.type = data_inhibitory\$my.cell.type), 
    )
    adata <- concat(c(data_excitatory\$T, data_inhibitory\$T))
    adata\$write_h5ad("RNA.h5ad")
    """
}

// Import metadata file
process import_metadata {
    publishDir "${params.baseOutputDir}"

    input:
    file(metadata_excitatory)
    file(metadata_inhibitory)

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


////////////////////////////////////////////////////////////////////////////////
// DNA data import and preprocessing
////////////////////////////////////////////////////////////////////////////////

// import sc-ATAC data
process import_atac {
    publishDir "${params.baseOutputDir}/ATAC/"
    cpus 2

    input:
    tuple val(datasetID), file(datasetFile)

    output:
    tuple val(datasetID), file("${datasetID}.h5ad")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    data = snap.pp.import_data(
        "${datasetFile}",
        "${file(gffFile)}",
        snap.genome.GRCm38,
        sorted_by_barcode=False,
        file = "${datasetID}.h5ad",
        min_num_fragments=1000,
        min_tsse=3,
    )
    """
}

// Plot basic QC metric
process qc_stat {
    publishDir "result/Figures/QC/"

    input:
    path(dataset)

    output:
    path("*.pdf")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    import seaborn as sns
    import pandas as pd
    from matplotlib import pyplot as plt

    data = snap.read_dataset("${dataset}", no_check=True, mode = 'r')
    samples = np.unique(data.obs['sample'])
    dfs = []
    for sample in samples: 
        select = data.obs['sample'] == sample 
        df = pd.concat([
            pd.DataFrame({
                "value": data.adatas.obs['tsse'][select],
                "metric": "TSSe",
            }),
            pd.DataFrame({
                "value": data.adatas.obs['n_fragment'][select],
                "metric": "#Fragment",
            }),
            pd.DataFrame({
                "value": data.adatas.obs['frac_dup'][select],
                "metric": "fraction duplicated",
            }),
        ])
        df['sample'] = sample
        dfs.append(df)
    data.close()
    df = pd.concat(dfs)

    #sns.set_theme(style="whitegrid")
    #sns.set(rc={'figure.figsize':(7.2,7.2)})
    g = sns.catplot(
        data=df, x='sample', y='value',
        row='metric', kind='violin', sharey=False,
    )
    g.set_xticklabels(rotation=90)
    plt.savefig("atac_qc.pdf")
    """
}

// Merge all sc-ATAC h5ad files
process merge_data {
    publishDir "${params.baseOutputDir}/ATAC/"

    input:
    val(dataset)
    file(metadata)

    output:
    file("merged/_dataset.h5ads")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import polars as pl
    import numpy as np
    # merge .h5ad files
    files = ${dataset}
    data = snap.create_dataset(
        [(files[i], files[i+1]) for i in range(0, len(files), 2)],
        storage="merged.h5ads",
    )

    # subset dataset
    metadata = pl.read_csv("${metadata}", sep='\t')
    metadata = metadata[metadata['my.cell.type'].to_numpy() != None]

    obs = data.obs[:]
    obs.insert_at_idx(0, pl.Series("ID", obs['sample'] + "+" + obs['Cell']))
    data.obs = obs

    join = obs.join(metadata, on = "ID")
    common_barcodes = set(join["ID"])
    data_subset = data.subset(
        [i in common_barcodes for i in obs["ID"]],
        out = "merged",
    )
    data.close()

    idx_map = dict(zip(
        list(data_subset.obs["ID"]), range(data_subset.shape[0])
    ))
    order = np.argsort(np.array([idx_map[i] for i in join['ID']]))
    new_obs = join[order, :]
    new_obs = new_obs[:, [s.null_count() == 0 for s in new_obs]]

    assert np.all(new_obs["ID"].to_numpy() == data_subset.obs["ID"])
    data_subset.obs = new_obs
    data_subset.close()
    """
}

/*
process call_peaks {
    input:
    val(dataset)

    output:
    val(dataset)

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import os
    #os.chdir("${dataset}")
    data = snap.read_dataset("${dataset}/_dataset.h5ads", no_check = True)
    snap.tl.call_peaks(data, group_by = "my.cell.type", q_value = 0.005)
    """
}

process make_cell_by_peak_mat {
    publishDir 'result'
    cache 'lenient'
    input:
    val(dataset)

    output:
    path("peak_matrix.h5ad")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap

    data = snap.read_dataset("${dataset}/_dataset.h5ads", no_check = True)
    snap.pp.make_peak_matrix(data, file = "peak_matrix.h5ad")
    """
}

process diff {
    publishDir 'result'
    cache 'lenient'
    input:
    val(dataset)

    output:
    path("peak_matrix.h5ad")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    data = snap.read_dataset("${dataset}/_dataset.h5ads", no_check = True)
    snap.pp.make_peak_matrix(data, file = "peak_matrix.h5ad")
    """
}

process export_counts {
    publishDir 'result'
    cache 'lenient'
    input:
    path(peak_matrix)

    output:
    path("accessibility.tsv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import pandas as pd
    import numpy as np 
    data = snap.read("${peak_matrix}", mode = "r")
    rpkm = snap.tl.aggregate_X(
        data,
        group_by = "my.cell.type",
        normalize = "RPKM",
        inplace = False
    )
    df = pd.DataFrame(rpkm)
    df.index = data.var_names
    df = np.log2(df + 1)
    df.to_csv("accessibility.tsv", sep='\t')
    """
}

workflow {
    prepare_rna_data()
    metadata = prepare_meta_data()

    make_cell_by_peak_mat(call_peaks(merged))
}
*/


////////////////////////////////////////////////////////////////////////////////
// Workflow definition
////////////////////////////////////////////////////////////////////////////////

workflow import_data {
    take: data
    /*
    data.rds_file_excitatory
    data.rds_file_inhibitory
    data.metadata_excitatory
    data.metadata_inhibitor
    data.atac_fragment_files
    */

    main:
        import_rna(data.rds_file_excitatory, data.rds_file_inhibitory)

        metadata = import_metadata(data.metadata_excitatory, data.metadata_inhibitory)

        atac_data = data.atac_fragment_files.map { file ->
            tuple(file.toString().split("/").reverse()[1], file)
        } | import_atac

        merge_data(
            atac_data.collect { a, b -> ["\"${a}\"", "\"${b}\""] },
            metadata
        ) | qc_stat

    emit:
        merged_atac_data = merge_data.out
}