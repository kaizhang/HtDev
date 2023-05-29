nextflow.enable.dsl=2

params.baseOutputDir = 'result/atac_peaks/'

process call_peaks {
    input:
    path(merged_data)

    output:
    path("peaks.tsv.gz")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    data = snap.read_dataset("${merged_data}", mode='r')
    peaks = snap.tl.call_peaks(data, groupby="my.cell.type", q_value=0.005, inplace=False)
    data.close()
    peaks.write_csv("peaks.tsv.gz", sep = '\t')
    """
}

process output_peaks {
    publishDir 'result/peaks/'
    input:
    path("peaks.tsv.gz")

    output:
    path("*.bed.gz")

    """
    #!/usr/bin/env python3
    import polars as pl
    import snapatac2 as snap
    import gzip
    df = pl.read_csv("peaks.tsv.gz", sep = '\t')
    for name in df.columns[1:]:
        with gzip.open(name + ".bed.gz", "wt") as f:
            f.write(
                '\\n'.join([x.replace("-", "\t").replace(":", "\t") for x in df['Peaks'][df[name]]])
            )
    """
}

process make_cell_by_peak_mat {
    publishDir 'result/peaks/'

    input:
    path(merged_data)
    path("peaks.tsv")

    output:
    path("peak_matrix.h5ad")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    with open("peaks.tsv", 'r') as f:
        f.readline()
        peaks = [line.strip().split('\t')[0] for line in f]
    data = snap.read_dataset("${merged_data}", mode='r')
    snap.pp.make_peak_matrix(data, use_rep=peaks, file="peak_matrix.h5ad").close()
    data.close()
    """
}

process umap_embedding {
    publishDir 'result/Figures/QC/'

    input:
    path("peak_matrix.h5ad")

    output:
    path("*.html")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    import gzip

    data = snap.read("peak_matrix.h5ad", backed=None)
    snap.tl.spectral(data, features = None)
    snap.tl.umap(data)
    data.obs['neuron_type'] = [x.split('-')[0] for x in data.obs['my.cell.type']]
    snap.pl.umap(data, color = 'neuron_type', out_file="ATAC_embedding_neuron_type.html")
    snap.pl.umap(data, color = 'age', out_file = "ATAC_embedding_age.html")
    """
}

////////////////////////////////////////////////////////////////////////////////
// Workflow definition
////////////////////////////////////////////////////////////////////////////////

workflow atac_call_peaks {
    take: data
    /*
    data.merged_dataset
    */

    main:
        peaks = call_peaks(data.merged_dataset)
        output_peaks(peaks)
        umap_embedding(
            make_cell_by_peak_mat(data.merged_dataset, peaks)
        )

    emit:
        peak_matrix = make_cell_by_peak_mat.out
}
