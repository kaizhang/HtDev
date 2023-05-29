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
    peak_mat = snap.read("peak_mat.h5ad", backed='r')
    snap.tl.aggregate_X(
        peak_mat,
        groupby='my.cell.type',
        normalize="RPKM",
        file="cell_type_accessibility.h5ad",
    ).close()
    peak_mat.close()
    """
}

process specific_peaks {
    input:
    path("peak_mat.h5ad")

    output:
    path("specific_peaks.h5ad")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    from scipy.stats import zscore
    import numpy as np
    import scipy.stats
    from statsmodels.stats.multitest import multipletests
    from scipy.sparse import csr_matrix

    def zscores_to_fdr(zs):
        return multipletests(scipy.stats.norm.sf(zs), alpha = 0.05, method = "fdr_bh")

    peak_mat = snap.read("peak_mat.h5ad", backed=None)

    mat = np.log2(peak_mat.X[:] + 1)
    z = zscore(mat, axis = 0)

    peaks = []
    for i in range(z.shape[0]):
        select = zscores_to_fdr(z[i, :])[0]
        peaks.append(select)

    csr = csr_matrix(np.vstack(peaks), dtype = np.float64)
    peak_mat.layers["specific_peaks"] = csr

    peak_mat.write("specific_peaks.h5ad", compression = "gzip")
    """
}

process plot_specific_peaks {
    publishDir "result/Figures/"

    input:
    path("peak_mat.h5ad")

    output:
    path("*.pdf")

    """
    #!/usr/bin/env Rscript
    library(anndata)
    library(ComplexHeatmap)
    library(viridis)
    data <- read_h5ad("peak_mat.h5ad")
    peaks <- data\$layers['specific_peaks']
    indices <- unique(peaks@j+1)
    data <- data[,indices]
    pdf("specific_peaks.pdf")
    draw(Heatmap(
        log2(t(data\$X) + 1),
        col=viridis(100),
        name = "log2(1+RPKM)",
        use_raster=T,
        cluster_rows=F,
        cluster_columns=F,
        clustering_method_columns="ward.D2",
        show_row_names=F,
        column_names_gp = gpar(fontsize=6),
        column_title = "cell types",
        row_title = "peaks",
    ))
    dev.off()
    """
}

process output_specific_peaks {
    input:
    path("data.h5ad")

    output:
    path("*.bed.gz")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import gzip
    data = snap.read("data.h5ad", backed=None)
    peaks = data.layers["specific_peaks"]
    for i, name in enumerate(data.obs_names):
        ps = data.var_names[peaks[i, :].indices]
        if len(ps) > 100:
            with gzip.open(name + ".bed.gz", "wt") as f:
                f.write("\\n".join([p.replace(":", "\t").replace("-", "\t") for p in ps]))
    """
}

process specific_peak_go {
    input:
    path(peak)

    output:
    path("*.tsv")

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
    write.table(tb, file=paste0(strsplit("${peak}"[[1]][1], ".bed.gz"), "_GO.tsv"), sep="\t", quote=F, row.names=F)
    """
}

process output_go {
    publishDir "${params.baseOutputDir}/specific_peaks/"

    input:
    path(files)

    output:
    path("*.tsv")

    """
    #!/usr/bin/env python
    import pandas as pd
    files = "${files}".split(" ")
    terms = set()
    filtered_files = []
    for fl in files:
        df = pd.read_csv(fl, sep="\t")
        df = df[(df.Binom_Fold_Enrichment >= 2) & (df.Binom_Adjp_BH <= 0.05)]
        df.sort_values("Binom_Adjp_BH", inplace=True)
        df = df.head(10)
        for term in df['name']:
            terms.add(term)
        if df.shape[0] > 0:
            filtered_files.append(fl)

    data = {}
    for fl in filtered_files:
        name = fl.split("/")[-1].split("_")[0]
        df = pd.read_csv(fl, sep="\t", index_col=1)
        data[name] = df.loc[terms]['Binom_Adjp_BH']
    pd.DataFrame(data).fillna(1).to_csv("GO_enrichment.tsv", sep="\t")
    """
}

process plot_go {
    publishDir "result/Figures/"
    input:
    path("GO_enrichment.tsv")

    output:
    path("*.pdf")

    """
    #!/usr/bin/env Rscript
    library("ggplot2")
    library("ComplexHeatmap")
    library("RColorBrewer")
    library("viridis")
    mat <- read.table("GO_enrichment.tsv", sep="\t", header=T, row.names=1)
    mat <- -log10(mat)
    mat[mat>30] <- 30
    pdf("restricted_peaks_go.pdf", width=8.0, height=20.0)
    draw(Heatmap(mat, 
        col=magma(100),
        cluster_rows=T, cluster_columns=T,
        show_row_dend=F,
        show_column_dend=F,
        show_row_names=T, row_names_side="left",
        show_column_names=T,
        row_names_gp = gpar(fontsize=5),
        column_names_gp = gpar(fontsize=5),
        heatmap_legend_param = list(
            title = "-log10(P-value)",
            title_position = "leftcenter-rot"
            #at = c(-2, 0, 2)
        )
    ))
    dev.off()
    """
}

process specific_peak_motif {
    input:
    path("data.h5ad")

    output:
    path("*.csv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import gzip
    data = snap.read("data.h5ad", backed=None)

    peaks = data.layers["specific_peaks"]
    peak_list = {}
    for i, name in enumerate(data.obs_names):
        ps = data.var_names[peaks[i, :].indices]
        if len(ps) > 100:
            peak_list[name] = ps
    result = snap.tl.motif_enrichment(
        motifs=snap.datasets.Meuleman_2020(),
        regions=peak_list,
        genome_fasta=snap.genome.mm10,
    )
    for k, df in result.items():
        df.sort('adjusted p-value').write_csv(k + ".csv")
    """
}

process output_motif {
    publishDir "${params.baseOutputDir}/specific_peaks/"
    input:
    path(files)

    output:
    path("*.tsv")

    """
    #!/usr/bin/env python
    import pandas as pd
    files = "${files}".split(" ")
    df = pd.read_csv(files[0])
    tf = {}
    for i, row in df.iterrows():
        tf[row['family']] = row['name']

    dfs = []
    for fl in files:
        name = fl.split("/")[-1].split(".")[0]
        df = pd.read_csv(fl)
        df = df[df['adjusted p-value'] <= 0.05]
        df = df.groupby("family").min().sort_values("adjusted p-value")[['adjusted p-value']]
        df.columns = [name]
        if df.shape[0] > 0:
            dfs.append(df.head(10))
    df = pd.concat(dfs, axis=1)
    df.index = [tf[i] for i in df.index]
    df = df[~df.index.duplicated(keep='first')]
    pd.DataFrame(df).fillna(1).to_csv("motif_enrichment.tsv", sep="\t")
    """
}

process plot_motif {
    publishDir "result/Figures/"
    input:
    path("motif_enrichment.tsv")

    output:
    path("*.pdf")

    """
    #!/usr/bin/env Rscript
    library("ggplot2")
    library("ComplexHeatmap")
    library("RColorBrewer")
    library("viridis")
    mat <- read.table("motif_enrichment.tsv", sep="\t", header=T, row.names=1)
    mat <- -log10(mat)
    mat[mat>30] <- 30
    pdf("restricted_peaks_motif.pdf", width=8.0, height=14.0)
    draw(Heatmap(mat, 
        col=magma(100),
        cluster_rows=T, cluster_columns=T,
        show_row_dend=F,
        show_column_dend=F,
        show_row_names=T, row_names_side="left",
        show_column_names=T,
        row_names_gp = gpar(fontsize=5),
        column_names_gp = gpar(fontsize=5),
        heatmap_legend_param = list(
            title = "-log10(P-value)",
            title_position = "leftcenter-rot"
            #at = c(-2, 0, 2)
        )
    ))
    dev.off()
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
        plot_specific_peaks(sp_peak)

        output_specific_peaks(sp_peak)
            | flatten | specific_peak_go | collect | output_go | plot_go

        specific_peak_motif(sp_peak) | output_motif | plot_motif
}
