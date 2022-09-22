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
    import snapatac2 as snap
    peak_mat = snap.read("${peak_mat}", mode = 'r')
    snap.tl.aggregate_X(
        peak_mat,
        group_by = peak_mat.obs['my.cell.type'] + ":" + peak_mat.obs['age'],
        normalize = "RPM",
        file = "peak_matrix_aggr.h5ad",
    ).close()
    peak_mat.close()

    gene_mat = snap.read("${gene_mat}", mode = 'r')
    snap.tl.aggregate_X(
        gene_mat,
        group_by = gene_mat.obs['my.cell.type'] + ":" + gene_mat.obs['age'],
        normalize = "RPM",
        file = "gene_matrix_aggr.h5ad",
    ).close()
    gene_mat.close()
    """
}

process init_network {
    input:
    tuple path(peak_mat), path(gene_mat)

    output:
    file("network.p")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import gzip
    from pickle import dump
    peak_mat = snap.read("${peak_mat}", mode = "r")
    gene_mat = snap.read("${gene_mat}", mode = "r")
    network = snap.tl.init_network_from_annotation(
        peak_mat.var_names,
        "${file(gffFile)}",
        upstream=500000,
        downstream=500000,
    )
    snap.tl.add_cor_scores(network, peak_mat, gene_mat)
    dump(network, gzip.open("network.p", "wb"))
    """
}

process perform_regression {
    input:
    tuple path("peak_mat.h5ad"), path("gene_mat.h5ad")
    path("input")

    output:
    path("network.p")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import gzip
    from pickle import dump, load
    network = load(gzip.open("input", "rb"))
    peak_mat = snap.read("peak_mat.h5ad", mode = "r")
    gene_mat = snap.read("gene_mat.h5ad", mode = "r")
    snap.tl.prune_network(
        network,
        edge_filter = lambda x: abs(x.correlation_score) > 0.05,
    )
    snap.tl.add_regr_scores(network, peak_mat, gene_mat, use_gpu=True)
    dump(network, gzip.open("network.p", "wb"))
    """
}

process find_motifs {
    input:
    path("input")

    output:
    path("network.p")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    import gzip
    from pickle import dump, load
    network = load(gzip.open("input", "rb"))
    cutoff = np.quantile([e.regr_score for e in network.edges()], 0.9)
    snap.tl.prune_network(
        network,
        edge_filter = lambda edge: edge.regr_score > cutoff and abs(edge.correlation_score) > 0.05,
    )
    snap.tl.add_tf_binding(
        network = network,
        motifs = snap.datasets.cis_bp(unique = True),
        genome_fasta = snap.genome.mm10,
    )
    dump(network, gzip.open("network.p", "wb"))
    """
}

process network_stat {
    publishDir "result/Figures/network/"
    input:
    path("input")

    output:
    path("*.pdf")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import gzip
    from pickle import dump, load
    network = load(gzip.open("input", "rb"))
    snap.pl.network_scores(network, "correlation_score", interactive = False, out_file = "correlation.pdf")
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
        matrix_file = aggregate_cells(data.peak_matrix, data.gene_matrix)
        network = perform_regression(
            matrix_file,
            init_network(matrix_file)
        )
        network_stat(network)
        find_motifs(network)
}
