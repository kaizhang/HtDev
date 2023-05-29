nextflow.enable.dsl=2

process list_regions {
    input:
    path("data.h5ad")

    output:
    path("*")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    data = snap.read("data.h5ad", backed='r')
    labels = np.unique(data.obs['brain.region'])
    for region in labels:
        with open(region, "w") as f:
            f.write("")
    """
}

process region_diff {
    input:
    tuple val(label), path("data.h5ad")

    output:
    path("*.tsv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    data = snap.read("data.h5ad", backed=None)
    region = "${label}".split("/")[-1]
    fg = np.where(data.obs['brain.region'] == region)[0]
    bg = np.where(data.obs['brain.region'] != region)[0]
    diff_peaks = snap.tl.diff_test(
        data,
        cell_group1=fg,
        cell_group2=bg,
        direction='positive',
    )
    diff_peaks.write_csv(region + ".tsv", sep="\t")
    """
}

process specific_peak_matrix {
    input:
    path(files)
    path("peak_mat.h5ad")

    output:
    path("*.h5ad")

    """
    #!/usr/bin/env python3
    import pandas as pd
    import anndata as ad
    import snapatac2 as snap
    files = ${files.collect {"'" + it + "'"}}
    peak_mat = snap.read("peak_mat.h5ad", backed=None)
    peak_mat = snap.tl.aggregate_X(peak_mat, groupby="brain.region", normalize="RPKM")
    adatas = []
    for f in files:
        name = f.split("/")[-1].split(".tsv")[0]
        df = pd.read_csv(f, sep="\t")
        idx = df["adjusted p-value"] <= 0.01
        features = df[idx].index
        data = peak_mat[:, features].copy().T
        data.obs['brain.region'] = name
        adatas.append(data)
    ad.concat(adatas).write_h5ad("specific_peaks.h5ad")
    """
}

process specific_peak_matrix2 {
    input:
    path("peak_mat.h5ad")

    output:
    path("*.h5ad")

    """
    #!/usr/bin/env python3
    import pandas as pd
    import anndata as ad
    import snapatac2 as snap

    peak_mat = snap.read("peak_mat.h5ad", backed=None)
    peaks = snap.tl.marker_regions(peak_mat, groupby="brain.region", pvalue=0.001)
    peak_mat = snap.tl.aggregate_X(peak_mat, groupby="brain.region", normalize="RPKM")

    adatas = []
    for name, p in peaks.items():
        data = peak_mat[:, p].copy().T
        data.obs['brain.region'] = name
        adatas.append(data)
    ad.concat(adatas).write_h5ad("specific_peaks.h5ad")
    """
}

process plot_specific_peak {
    publishDir "result/Figures/"

    input:
    path("peak_mat.h5ad")

    output:
    path("*.pdf")

    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    import snapatac2 as snap
    import matplotlib.pylab as plt
    import PyComplexHeatmap
    from PyComplexHeatmap import *
    import glasbey

    data = snap.read("peak_mat.h5ad", backed=None)
    df = pd.DataFrame(np.log(data.X + 1), index=data.obs_names, columns=data.var_names)

    plt.figure(figsize=(5, 7))
    regions = np.unique(data.obs['brain.region'])
    colors = glasbey.create_palette(palette_size=regions.size)

    row_ha_left = HeatmapAnnotation(
        label=anno_label(data.obs['brain.region'], merge=True, extend=True, colors=colors),
        Group=anno_simple(
            data.obs['brain.region'],
            colors=colors,
            legend=False,
        ),
        axis=0,verbose=0,orientation='left'
    )
    cm = ClusterMapPlotter(
        data=df, 
        left_annotation=row_ha_left,
        rasterized=True,
        row_split=data.obs['brain.region'],
        row_split_gap=5,
        row_cluster=False,
        col_cluster=False,
        show_colnames=True,
        #cmap='viridis',
        legend_gap=5,legend_width=5,legend_hpad=2,legend_vpad=5,
    )
    plt.savefig('region_peaks.pdf', bbox_inches='tight')
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
    import numpy as np
    import gzip
    data = snap.read("data.h5ad", backed=None)
    regions = np.unique(data.obs['brain.region'])
    for name in regions:
        with gzip.open(name + ".bed.gz", "wt") as f:
            f.write("\\n".join([p.replace(":", "\t").replace("-", "\t") for p in data[data.obs['brain.region'] == name].obs_names]))
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
        df = df.head(5)
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
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    import matplotlib.pylab as plt
    import PyComplexHeatmap
    from PyComplexHeatmap import *
    import glasbey

    df = pd.read_csv("GO_enrichment.tsv", sep="\t", index_col=0)
    df = -np.log10(df)
    df = df.clip(upper=30)
    df = df.rename(columns={"AvPe.MnPO": "AvPe.MnPO.inhibitory", "AvPE.MnPO": "AvPe.MnPO.excitatory"})

    plt.figure(figsize=(5, 20))
    cm = ClusterMapPlotter(
        data=df, 
        row_cluster=True,
        col_cluster=True,
        show_colnames=True,
        show_rownames=True,
        cmap='viridis',
        legend_gap=5,legend_width=5,legend_hpad=2,legend_vpad=5,
    )
    plt.savefig('region_peaks_go.pdf', bbox_inches='tight')
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
    import numpy as np
    data = snap.read("data.h5ad", backed=None)
    regions = np.unique(data.obs['brain.region'])
    peak_list = {}
    for name in regions:
        peak_list[name] = [p for p in data[data.obs['brain.region'] == name].obs_names]
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
        name = fl.split("/")[-1].split(".csv")[0]
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
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    import matplotlib.pylab as plt
    import PyComplexHeatmap
    from PyComplexHeatmap import *
    import glasbey

    df = pd.read_csv("motif_enrichment.tsv", sep="\t", index_col=0)
    df = -np.log10(df)
    df = df.clip(upper=30)
    df = df.rename(columns={"AvPe.MnPO": "AvPe.MnPO.inhibitory", "AvPE.MnPO": "AvPe.MnPO.excitatory"})

    plt.figure(figsize=(4, 18))
    cm = ClusterMapPlotter(
        data=df, 
        row_cluster=True,
        col_cluster=True,
        show_colnames=True,
        show_rownames=True,
        cmap='viridis',
        legend_gap=5,legend_width=5,legend_hpad=2,legend_vpad=5,
        xticklabels_kws={'labelsize':6},
        yticklabels_kws={'labelsize':6},
    )
    plt.savefig('region_peaks_motif.pdf', bbox_inches='tight')
    """
}

////////////////////////////////////////////////////////////////////////////////
// Workflow definition
////////////////////////////////////////////////////////////////////////////////

workflow regional_analysis {
    take: data
    /*
    data.peak_matrix
    */

    main:
        peaks = list_regions(data.peak_matrix) | flatten | combine(data.peak_matrix)
            | region_diff

        specific_peak_matrix(peaks | collect, data.peak_matrix)
            //| plot_specific_peak

        specific_peaks = specific_peak_matrix2(data.peak_matrix)
        specific_peaks | plot_specific_peak
        specific_peaks | output_specific_peaks
            | flatten | specific_peak_go | collect | output_go | plot_go

        specific_peaks | specific_peak_motif | output_motif | plot_motif
}
