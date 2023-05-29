nextflow.enable.dsl=2

process download_GO {
    output:
    tuple path("go-basic.obo"), path("gene2go")

    """
    #!/usr/bin/env python3
    from goatools.base import download_go_basic_obo
    from goatools.base import download_ncbi_associations
    obo_fname = "go-basic.obo"
    gaf_file = "gene2go"
    download_go_basic_obo(obo_fname)
    download_ncbi_associations(gaf_file)
    """
}

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

process list_celltypes {
    input:
    path("data.h5ad")

    output:
    path("*")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    data = snap.read("data.h5ad", backed='r')
    labels = np.unique(data.obs['my.cell.type'])
    for region in labels:
        with open(region, "w") as f:
            f.write("")
    """
}

process split_by_region {
    input:
    tuple val(region), path("atac.h5ad"), path("rna.h5ad")

    output:
    tuple val(region), path("atac_out.h5ad"), path("rna_out.h5ad")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    atac = snap.read("atac.h5ad", backed='r')
    atac.subset(atac.obs['brain.region'] == "${region}", out='atac_out.h5ad')
    rna = snap.read("rna.h5ad", backed='r')
    rna.subset(rna.obs['brain.region'] == "${region}", out="rna_out.h5ad")
    """
}

process pls {
    input:
    tuple val(region), path("atac.h5ad"), path("rna.h5ad")

    output:
    tuple val(region), path("pls.h5ad")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import anndata as ad
    from scipy.sparse import csr_matrix
    from scipy.stats import zscore
    import numpy as np
    from sklearn.cross_decomposition import PLSCanonical

    def aggregate_counts(adata):
        data = snap.tl.aggregate_X(adata, groupby=[str(x) for x in adata.obs['age']], normalize='RPM')
        data.X = np.log1p(data.X)
        mask = [np.any(data.X[:, i] > 1) for i in range(data.X.shape[1])]
        data = data[:, mask].copy()
        data.X = zscore(data.X, axis=0)
        mask = [np.any(np.abs(data.X[:, i]) > 1.45) for i in range(data.X.shape[1])]
        data = data[:, mask].copy().T
        return data

    atac = snap.read("atac.h5ad", backed=None)
    atac = aggregate_counts(atac)
    rna = snap.read("rna.h5ad", backed=None)
    rna = aggregate_counts(rna)

    data = snap.AnnData(X=np.vstack([atac.X, rna.X]), filename="pls.h5ad")
    data.var_names = atac.var_names
    data.obs_names = list(atac.obs_names) + list(rna.obs_names)
    data.obs['modality'] = ['atac'] * atac.X.shape[0] + ['rna'] * rna.X.shape[0]
    data.close()
    """
}

process kmeans {
    input:
    tuple val(region), path("data.h5ad")

    output:
    tuple val(region), path("out.h5ad")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    from sklearn.cluster import KMeans
    import numpy as np
    import pandas as pd
    import matplotlib.pylab as plt
    import PyComplexHeatmap
    from PyComplexHeatmap import *
    data = snap.read("data.h5ad", backed=None)

    def find_elbow_point(X):
        inertias = []
        for k in range(4, 20):
            kmeans = KMeans(n_clusters=k, init='k-means++', max_iter=300, n_init=10, random_state=0)
            kmeans.fit(X)
            inertias.append(kmeans.inertia_)
    
        # Calculate the first derivative of the inertia curve
        dy = np.diff(inertias)
        dx = np.diff(range(4, 20))
        gradient = dy / dx

        # Find the elbow point
        elbow_point = np.argmax(gradient) + 2

        return elbow_point

    n = find_elbow_point(data.X)

    kmeans = KMeans(n_clusters=n, init='k-means++', max_iter=300, n_init=10, random_state=0).fit(data.X)
    cluster_labels = [str(x) for x in kmeans.labels_]

    data.obs['clusters'] = cluster_labels
    data.write("out.h5ad", compression='gzip')
    """
}

process filter_peaks {
    input:
    tuple val(region), path("data.h5ad")

    output:
    tuple val(region), path("filtered.h5ad")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np

    def filter_peaks(peaks, genes):
        gene_set = {g.upper() for g in genes}
        network = snap.tl.init_network_from_annotation(
            peaks,
            snap.genome.mm10,
            upstream=500000,
            downstream=500000,
        )
        network = snap.tl.prune_network(
            network,
            node_filter=lambda x: x.type == 'region' or x.id in gene_set)
        regions = {node.id for node in network.nodes() if node.type == 'region'}
        peaks = {p for p in peaks if p in regions}
        peaks.update(genes)
        return peaks

    data = snap.read("data.h5ad", backed='r')
    clusters = data.obs['clusters'].unique()
    features = set()
    for cluster in clusters:
        peaks = (data.obs['clusters'].to_numpy() == cluster) & (data.obs['modality'].to_numpy() == 'atac')
        peaks = np.array(data.obs_names)[peaks]
        genes = (data.obs['clusters'].to_numpy() == cluster) & (data.obs['modality'].to_numpy() == 'rna')
        genes = np.array(data.obs_names)[genes]
        features.update(filter_peaks(peaks, genes))
    mask = [x in features for x in data.obs_names]
    data.subset(mask, out="filtered.h5ad")
    """
}

process motif_enrichment {
    publishDir "result/by_stage/motif/${region}/"
    input:
    tuple val(region), path("data.h5ad")
    path("all.h5ad")

    output:
    path("*.csv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    import gzip
    data = snap.read("data.h5ad", backed=None)
    # random sample
    all_peaks = np.array(snap.read('all.h5ad', backed=None).var_names)
    background = all_peaks[np.random.choice(all_peaks.size, 10000, replace=False)]

    clusters = data.obs['clusters'].unique()
    peak_list = {}
    for cluster in clusters:
        peaks = (data.obs['clusters'].to_numpy() == cluster) & (data.obs['modality'].to_numpy() == 'atac')
        peak_list[cluster] = np.array(data.obs_names)[peaks]
    result = snap.tl.motif_enrichment(
        motifs=snap.datasets.Meuleman_2020(),
        regions=peak_list,
        genome_fasta=snap.genome.mm10,
        background=background,
    )
    for k, df in result.items():
        df.sort('adjusted p-value').write_csv(k + ".csv")
    """
}

process go_enrichment {
    publishDir "result/by_stage/GO/${region}/"
    input:
    tuple val(region), path("data.h5ad")
    tuple path(obo_fname), path(gaf_file)

    output:
    path("*.tsv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    import gzip
    import polars as pl
    from goatools.obo_parser import GODag
    from goatools.associations import read_ncbi_gene2go
    from goatools.test_data.genes_NCBI_10090_ProteinCoding import GENEID2NT as geneid2nt_mouse
    from goatools.go_enrichment import GOEnrichmentStudy

    geneSymbol2id = {v.Symbol: k for k, v in geneid2nt_mouse.items()}

    obodag = GODag("${obo_fname}")
    geneid2gos_mouse = read_ncbi_gene2go("${gaf_file}", taxids=[10090])

    data = snap.read("data.h5ad", backed=None)

    clusters = data.obs['clusters'].unique()
    for cluster in clusters:
        genes = (data.obs['clusters'].to_numpy() == cluster) & (data.obs['modality'].to_numpy() == 'rna')
        genes = np.array(data.obs_names)[genes]
        gene_list = []
        for g in genes:
            if g in geneSymbol2id:
                gene_list.append(geneSymbol2id[g])
        # Run GO term enrichment analysis
        goeaobj = GOEnrichmentStudy(
            geneid2nt_mouse.keys(),  # All protein-coding genes
            geneid2gos_mouse,  # geneid/GO associations
            obodag,
            propagate_counts=False,
            alpha=0.05,  # FDR cutoff
            methods=['fdr_bh'],  # FDR correction method
        )
        goea_results_all = goeaobj.run_study(gene_list)
        # Filter significant GO terms
        goea_results_sig = [result for result in goea_results_all if result.p_fdr_bh < 0.05]
        # Convert significant GO terms to pandas DataFrame
        goea_df = pl.DataFrame(
            [
                {
                    "GO_term": result.GO,
                    "name": result.name,
                    "enrichment": result.enrichment,
                    "p_fdr_bh": result.p_fdr_bh,
                    "p_uncorrected": result.p_uncorrected,
                }
                for result in goea_results_sig
            ]
        )
        if goea_df.shape[0] > 1:
            geoa_df = goea_df.sort('p_fdr_bh')
        goea_df.write_csv(cluster + ".tsv", sep='\t')
    """
}




process plot_kmeans {
    publishDir "result/by_stage/heatmap/"

    input:
    tuple val(region), path("data.h5ad")

    output:
    path("*.pdf")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    from sklearn.cluster import KMeans
    import numpy as np
    import pandas as pd
    import matplotlib.pylab as plt
    import PyComplexHeatmap
    from PyComplexHeatmap import *

    data = snap.read("data.h5ad", backed=None)

    meta = pd.DataFrame({
        'clusters': [str(x) for x in data.obs['clusters']],
        'modality': [str(x) for x in data.obs['modality']],
    }, index=data.obs_names)
    df = pd.DataFrame(data.X[:], columns=data.var_names, index=data.obs_names)

    plt.figure(figsize=(2, 4))
    row_ha = HeatmapAnnotation(
        label=anno_label(meta.clusters, merge=True),
        clusters=anno_simple(meta.clusters, rasterized=True),
        modality=anno_simple(meta.modality, rasterized=True),
        axis=0,
    )
    cm = ClusterMapPlotter(
        left_annotation=row_ha,
        data=df, 
        row_split=meta,
        row_split_gap=0.0,
        row_cluster=False,
        col_cluster=False,
        show_colnames=True,
        show_rownames=False,
        cmap='viridis',
        rasterized=True,
        legend_gap=5,legend_width=5,legend_hpad=2,legend_vpad=5,
        xticklabels_kws={'labelsize':6},
        yticklabels_kws={'labelsize':6},
    )
    plt.savefig("${region}.pdf", bbox_inches='tight')
    """
}

process split_by_cell_type {
    input:
    tuple val(cell), path("atac.h5ad"), path("rna.h5ad")

    output:
    tuple val(cell), path("atac_out.h5ad"), path("rna_out.h5ad")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    atac = snap.read("atac.h5ad", backed='r')
    atac.subset(atac.obs['my.cell.type'] == "${cell}", out='atac_out.h5ad')
    rna = snap.read("rna.h5ad", backed='r')
    rna.subset(rna.obs['my.cell.type'] == "${cell}", out="rna_out.h5ad")
    """
}



process k_means2 {
    publishDir "result/test/"
    input:
    tuple val(region), path("data.h5ad")

    output:
    path("*.pdf")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    from sklearn.cluster import KMeans
    import numpy as np
    import pandas as pd
    import matplotlib.pylab as plt
    import PyComplexHeatmap
    from PyComplexHeatmap import *
    data = snap.read("data.h5ad", backed=None)
    data = snap.tl.aggregate_X(data, groupby=[str(x) for x in data.obs['age']], normalize='RPKM')
    data.X = np.log1p(data.X)
    mask = [np.any(data.X[:, i] > 1) for i in range(data.X.shape[1])]
    data = data[:, mask].copy().T

    kmeans = KMeans(n_clusters=10, random_state=0).fit(data.X)
    cluster_labels = kmeans.labels_
    reorder = np.argsort(cluster_labels)

    df = pd.DataFrame(data.X[reorder, :], columns=data.var_names, index=data.obs_names[reorder])

    plt.figure(figsize=(2, 4))
    cm = ClusterMapPlotter(
        data=df, 
        row_cluster=False,
        col_cluster=False,
        show_colnames=True,
        show_rownames=False,
        cmap='viridis',
        rasterized=True,
        legend_gap=5,legend_width=5,legend_hpad=2,legend_vpad=5,
        xticklabels_kws={'labelsize':6},
        yticklabels_kws={'labelsize':6},
    )
    plt.savefig("${region}.pdf", bbox_inches='tight')
 
    """
}



////////////////////////////////////////////////////////////////////////////////
// Workflow definition
////////////////////////////////////////////////////////////////////////////////

workflow stage_analysis {
    take: data
    /*
    data.peak_matrix
    data.gene_matrix
    */

    main:
        go = download_GO()

        list_regions(data.peak_matrix) | flatten
            | map { it.toString().split("/").last() }
            | combine(data.peak_matrix) | combine (data.gene_matrix)
            | split_by_region


        peak_cluster = list_celltypes(data.peak_matrix) | flatten
            | map { it.toString().split("/").last() }
            | combine(data.peak_matrix) | combine (data.gene_matrix)
            | split_by_cell_type
            | pls
            | kmeans
            | filter_peaks

        go_enrichment(peak_cluster, go)
        motif_enrichment(peak_cluster, data.peak_matrix)
        peak_cluster | plot_kmeans

    /*
        Channel.fromList(["i-B1", "i-M1", "i-M2", "e-A1", "e-P2"])
            | combine(data.peak_matrix)
    */
}