nextflow.enable.dsl=2

params.baseOutputDir = 'result/differential_analysis/'

process cell_type_accessibility {
    publishDir "${params.baseOutputDir}"

    input:
    path("peak_mat.h5ad")

    output:
    path("cell_type_accessibility.h5ad")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    peak_mat = snap.read("peak_mat.h5ad", mode = 'r')
    snap.tl.aggregate_X(
        peak_mat,
        group_by = 'my.cell.type',
        normalize = "RPKM",
        file = "cell_type_accessibility.h5ad",
    ).close()
    peak_mat.close()
    """
}

process specific_peaks {
    input:
    path("peak_mat.h5ad")

    output:
    path("peaks.p")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    from scipy.stats import zscore
    import numpy as np
    import scipy.stats
    from statsmodels.stats.multitest import multipletests
    import gzip
    from pickle import dump

    def zscores_to_fdr(zs):
        return multipletests(scipy.stats.norm.sf(zs), alpha = 0.05, method = "fdr_bh")

    peak_mat = snap.read("peak_mat.h5ad", mode = 'r')
    cell_names = np.array(peak_mat.obs_names)
    peak_names = np.array(peak_mat.var_names)

    mat = np.log2(peak_mat.X[:] + 1)
    z = zscore(mat, axis = 0)

    peaks = {}
    for i in range(z.shape[0]):
        select = zscores_to_fdr(z[i, :])[0]
        if np.where(select)[0].size >= 100:
            peaks[cell_names[i]] = peak_names[select]

    dump(peaks, gzip.open("peaks.p", "wb"))
    peak_mat.close()
    """
}

process plot_specific_peaks {
    publishDir "${params.baseOutputDir}"

    input:
    path("peak_mat.h5ad")
    path("peaks.p")

    output:
    path("*.pdf")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import pandas as pd
    import numpy as np
    import seaborn as sns
    from scipy.cluster import hierarchy
    from matplotlib import pyplot as plt
    import gzip
    from pickle import load

    peaks = load(gzip.open("peaks.p", "rb"))
    peak_mat = snap.read("peak_mat.h5ad", mode = 'r')

    mat = np.log2(peak_mat.X[:] + 1)

    vals = []
    for ps in peaks.values():
        j = peak_mat.var_ix(list(ps))
        vals.append(mat[:, j].mean(axis = 1))

    vals = pd.DataFrame(
        np.array(vals)[:, peak_mat.obs_ix(list(peaks.keys()))],
        index = peaks.keys(),
        columns = peaks.keys(),
    )

    Z = hierarchy.ward(vals)
    ordering = hierarchy.leaves_list(hierarchy.optimal_leaf_ordering(Z, vals))
    sns.heatmap(vals.iloc[ordering, ordering])
    plt.savefig("specific_peaks.pdf")

    peak_mat.close()
    """
}

process output_specific_peaks {
    publishDir "${params.baseOutputDir}/specific_peaks/"

    input:
    path("peaks.p")

    output:
    path("*.bed.gz")

    """
    #!/usr/bin/env python3
    import gzip
    from pickle import load
    peaks = load(gzip.open("peaks.p", "rb"))
    for k, ps in peaks.items():
        with gzip.open(k + ".bed.gz", "wt") as f:
            f.write("\\n".join([p.replace(":", "\t").replace("-", "\t") for p in ps]))
    """
}

process specific_peak_go {
    publishDir "${params.baseOutputDir}/specific_peaks/"

    input:
    path(peak)

    output:
    file("*.tsv.gz")

    """
    #!/usr/bin/env Rscript
    library("rGREAT")

    bed <- as.data.frame(read.table("${peak}", header = F, sep="\t", stringsAsFactors=FALSE, quote=""))
    job <- submitGreatJob(
        bed,
        species = "mm10",
        version = "4.0",
        request_interval = 10
    )
    tb <- getEnrichmentTables(job, ontology = c("GO Biological Process"))
    tb <- tb\$`GO Biological Process`
    select <- tb\$Binom_Fold_Enrichment >= 2 & tb\$Binom_Adjp_BH < 0.01
    tb <- tb[select, ]
    tb <- tb[order(tb\$Binom_Adjp_BH), ][1:50, ]
    write.table(tb, file=paste0(strsplit("${peak}"[[1]][1], ".bed.gz"), "_GO.tsv.gz"), sep="\t", quote=F, row.names=F)
    """
}

process specific_peak_motif {
    publishDir "${params.baseOutputDir}/specific_peaks/motifs/"

    input:
    path("peaks.p")

    output:
    file("*.csv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import gzip
    from pickle import load
    peaks = load(gzip.open("peaks.p", "rb"))
    result = snap.tl.motif_enrichment(
        motifs = snap.datasets.cis_bp(unique = True),
        regions = peaks,
        genome_fasta = snap.genome.mm10,
    )
    for k, df in result.items():
        df.sort('adjusted p-value').write_csv(k + ".csv")
    """
}


////////////////////////////////////////////////////////////////////////////////
// Workflow definition
////////////////////////////////////////////////////////////////////////////////

workflow differential_analysis {
    take: data
    /*
    data.peak_matrix
    */

    main:
        peak_mat = cell_type_accessibility(data.peak_matrix) 
        sp_peak = specific_peaks(peak_mat)
        plot_specific_peaks(peak_mat, sp_peak)

        output_specific_peaks(sp_peak)
            | flatten | specific_peak_go

        specific_peak_motif(sp_peak)
}
