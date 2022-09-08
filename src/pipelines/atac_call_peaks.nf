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
    data = snap.read_dataset("${merged_data}", no_check = True, mode='r')
    peaks = snap.tl.call_peaks(data, group_by = "my.cell.type", q_value = 0.005, inplace=False)
    data.close()
    peaks.write_csv("peaks.tsv.gz", sep = '\t')
    """
}

process make_cell_by_peak_mat {
    publishDir 'result/peaks/'

    input:
    path(merged_data)
    path(peak_file)

    output:
    path("peak_matrix.h5ad")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import gzip
    with gzip.open('${peak_file}', 'rt') as f:
        f.readline()
        peaks = [line.strip().split('\t')[0] for line in f]
    data = snap.read_dataset("${merged_data}", mode = 'r', no_check = True)
    snap.pp.make_peak_matrix(data, use_rep = peaks, file = "peak_matrix.h5ad").close()
    data.close()
    """
}

/*
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
*/


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
        make_cell_by_peak_mat(data.merged_dataset, peaks)

    emit:
        peak_matrix = make_cell_by_peak_mat.out
}