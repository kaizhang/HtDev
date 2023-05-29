nextflow.enable.dsl=2

params.baseOutputDir = 'result/regulatory_network/'

process aggregate_cells {
    input:
    path(peak_mat)
    path(gene_mat)

    output:
    tuple path("peak_matrix_aggr.h5ad"), path("gene_matrix_aggr.h5ad")

    """
    #!/usr/bin/env python3
    import scanpy as sc
    import snapatac2 as snap

    #adata = sc.read("${gene_mat}")
    #sc.pp.filter_genes(adata, min_cells=3)
    #sc.pp.normalize_total(adata, target_sum=1e4)
    #sc.pp.log1p(adata)
    #sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    #adata = adata[:, adata.var.highly_variable]
    #sc.pp.scale(adata, max_value=10)
    #sc.tl.pca(adata)

    #pseudo_cells = snap.tl.aggregate_cells(adata.obsm['X_pca'])


    gene_mat = snap.read("${gene_mat}", backed='r')
    pseudo_cells = gene_mat.obs['my.cell.type'] + ":" + gene_mat.obs['age']
    gene_mat_aggr = snap.tl.aggregate_X(
        gene_mat,
        groupby=pseudo_cells,
        normalize="RPM",
        file="gene_matrix_aggr.h5ad",
    )
    gene_mat_aggr.subset(
        obs_indices=[x != '-1' for x in gene_mat_aggr.obs_names],
        var_indices=gene_mat_aggr.X[:].sum(axis=0) != 0,
    )
    #mapping = dict(list(zip(gene_mat.obs_names, pseudo_cells)))
    gene_mat_aggr.close()
    gene_mat.close()

    peak_mat = snap.read("${peak_mat}", backed='r')
    #pseudo_cells = [mapping[x] for x in peak_mat.obs_names]
    pseudo_cells = peak_mat.obs['my.cell.type'] + ":" + peak_mat.obs['age']
    peak_mat_aggr = snap.tl.aggregate_X(
        peak_mat,
        groupby=pseudo_cells,
        normalize="RPM",
        file="peak_matrix_aggr.h5ad",
    )
    peak_mat_aggr.subset(obs_indices=[x != '-1' for x in peak_mat_aggr.obs_names])
    peak_mat_aggr.close()
    peak_mat.close()
    """
}

process init_network {
    input:
    tuple path(peak_mat), path(gene_mat)

    output:
    tuple path("network.p"), path("pmat.h5ad"), path("gmat.h5ad")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    import gzip
    from pickle import dump
    peak_mat = snap.read("${peak_mat}", backed=False)
    gene_mat = snap.read("${gene_mat}", backed=False)

    peak_mat.X = np.log(peak_mat.X[:] + 1)
    gene_mat.X = np.log(gene_mat.X[:] + 1)

    idx = np.apply_along_axis(lambda x: np.any(x > 1.5), 0, gene_mat.X[:])
    gene_mat = gene_mat[:, idx]

    idx = np.apply_along_axis(lambda x: np.any(x > 1.5), 0, peak_mat.X[:])
    peak_mat = peak_mat[:, idx]

    network = snap.tl.init_network_from_annotation(
        peak_mat.var_names,
        snap.genome.mm10,
        upstream=500000,
        downstream=500000,
    )
    snap.tl.add_cor_scores(network, peak_mat=peak_mat, gene_mat=gene_mat)
    dump(network, gzip.open("network.p", "wb"))
    peak_mat.write("pmat.h5ad")
    gene_mat.write("gmat.h5ad")
    """
}

process perform_regression {
    input:
    tuple path("input.p"), path("peak_mat.h5ad"), path("gene_mat.h5ad")

    output:
    tuple path("network.p"), path("peak_mat.h5ad"), path("gene_mat.h5ad")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import gzip
    from pickle import dump, load
    network = load(gzip.open("input.p", "rb"))
    peak_mat = snap.read("peak_mat.h5ad", backed="r")
    gene_mat = snap.read("gene_mat.h5ad", backed="r")
    snap.tl.add_regr_scores(network, peak_mat=peak_mat, gene_mat=gene_mat, method='gb_tree', use_gpu=True)
    dump(network, gzip.open("network.p", "wb"))
    """
}

process find_motifs {
    input:
    tuple path("input.p"), path("peak_mat.h5ad"), path("gene_mat.h5ad")

    output:
    tuple path("network.p"), path("peak_mat.h5ad"), path("gene_mat.h5ad")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    import gzip
    from pickle import dump, load
    network = load(gzip.open("input.p", "rb"))
    cutoff = np.quantile([e.regr_score for e in network.edges() if e.regr_score is not None], 0.9)
    network = snap.tl.prune_network(
        network,
        edge_filter = lambda fr, to, edge: edge.regr_score is not None and edge.regr_score >= cutoff,
    )
    snap.tl.add_tf_binding(
        network=network,
        motifs=snap.datasets.cis_bp(),
        genome_fasta=snap.genome.mm10,
    )
    dump(network, gzip.open("network.p", "wb"))
    """
}

process make_cell_type_network {
    input:
    path("peak_mat.h5ad")
    path("gene_mat.h5ad")
    tuple path("input.p"), path("peak_mat_aggr.h5ad"), path("gene_mat_aggr.h5ad")

    output:
    path("*")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    import gzip
    from pickle import dump, load
    def get_vars(x, mat):
        mask = np.ravel(mat[x, :].X) >= 1.5
        return set(np.array(mat.var_names)[mask])

    network = load(gzip.open("input.p", "rb"))
    for nd in network.nodes():
        if nd.type == "motif":
            nd.id = nd.id.capitalize()

    gene_mat = snap.read("gene_mat.h5ad", backed="r")
    gene_mat = snap.tl.aggregate_X(
        gene_mat,
        groupby="my.cell.type",
        normalize="RPM",
    )
    gene_mat.X = np.log(gene_mat.X[:] + 1)

    peak_mat = snap.read("peak_mat.h5ad", backed="r")
    peak_mat = snap.tl.aggregate_X(
        peak_mat,
        groupby="my.cell.type",
        normalize="RPM",
    )
    peak_mat.X = np.log(peak_mat.X[:] + 1)

    cell_type_nets = {}
    for cell_type in peak_mat.obs_names:
        regions = get_vars(cell_type, peak_mat)
        genes = get_vars(cell_type, gene_mat)
        net = snap.tl.prune_network(
            network,
            node_filter=lambda x: (x.type == "region" and x.id in regions) or (x.type != "region" and x.id in genes),
        )
        dump(net, gzip.open(cell_type, "wb"))
    """
}

process make_gene_network {
    input:
    tuple path(file), path("input.p"), path("peak_mat.h5ad"), path("gene_mat.h5ad")

    output:
    path("*.p")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    import gzip
    from pickle import dump, load
    from pathlib import Path

    file_path = Path("${file}")
    network = load(gzip.open(file_path, "rb"))
    cell_type = file_path.stem + ".p"

    peak_mat = snap.read("peak_mat.h5ad", backed="r")
    gene_mat = snap.read("gene_mat.h5ad", backed="r")

    network = snap.tl.link_tf_to_gene(network)
    snap.tl.add_regr_scores(
        network,
        peak_mat=peak_mat,
        gene_mat=gene_mat,
        alpha=0.1,
        scale_X=True,
        scale_Y=True,
    )
    network = snap.tl.prune_network(
        network,
        edge_filter=lambda fr, to, x: x.regr_score is not None and abs(x.regr_score) > 0 ,
    )
    dump(network, gzip.open(cell_type, "wb"))
    """
}

process pagerank {
    input:
    path(file)

    output:
    path("*.txt")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    import gzip
    from pickle import dump, load
    from pathlib import Path

    file_path = Path("${file}")
    network = load(gzip.open(file_path, "rb"))
    cell_type = file_path.stem

    result = snap.tools._network.pagerank(network)
    with open(cell_type + ".txt", 'w') as f:
        f.write('\\n'.join([k + '\t' + str(v) for k, v in result]))
    """
}

process merge_pagerank {
    publishDir "result/Figures/network/"
    input:
    val(files)

    output:
    tuple path("*.tsv"), path("*.html")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    from pathlib import Path
    import numpy as np
    import polars as pl
    from scipy.stats import zscore

    scores = None
    for file in ${files}:
        cell_type = Path(file).stem
        tfs = []
        vals = []
        with open(file, 'r') as f:
            for line in f:
                items = line.strip().split('\t')
                tfs.append(items[0])
                vals.append(float(items[1]))
        df = pl.DataFrame({"TF": tfs, cell_type: vals})
        if scores is None:
            scores = df
        else:
            scores = scores.join(df, on="TF", how="outer")

    scores = scores.fill_null(0.00001)
    scores = scores.to_pandas()
    scores = scores.set_index('TF')
    scores = np.log(scores)
    scores = scores.apply(zscore, axis=1)
    scores.to_csv("pagerank.tsv", sep='\t')

    max_z = [np.amax(row[1]) for row in scores.iterrows()]
    scores = scores.iloc[[abs(x) >= 3.1 for x in max_z], :]
    fig = snap.pl.heatmap(
        np.clip(scores.to_numpy(), -5, 5),
        row_names=scores.index,
        column_names=scores.columns,
        colorscale='RdBu_r',
    )
    snap.pl._base.render_plot(
        fig,
        out_file="pagerank.html",
        show=False,
        height=1400,
        width=800,
    )
    """
}

process network_stat {
    publishDir "result/Figures/network/"
    input:
    tuple path("input.p"), path("peak_mat.h5ad"), path("gene_mat.h5ad")

    output:
    path("*.pdf")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import gzip
    from pickle import dump, load
    network = load(gzip.open("input.p", "rb"))
    snap.pl.network_scores(network, "cor_score", interactive = False, out_file = "cor_score.pdf")
    snap.pl.network_scores(network, "regr_score", interactive = False, out_file = "regr_score.pdf")
    """
}


////////////////////////////////////////////////////////////////////////////////
// Workflow definition
////////////////////////////////////////////////////////////////////////////////

workflow regulatory_network {
    take: data
    /*
    data.peak_matrix
    data.gene_matrix
    */

    main:
        network = aggregate_cells(data.peak_matrix, data.gene_matrix) |
            init_network | perform_regression
        network_stat(network)
        make_cell_type_network(
            data.peak_matrix,
            data.gene_matrix,
            find_motifs(network),
        ) | flatten | combine(network) | make_gene_network | pagerank |
            map { "\"" + it + "\""} | collect | merge_pagerank
}