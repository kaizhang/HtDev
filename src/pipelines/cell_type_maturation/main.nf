nextflow.enable.dsl=2

process list_celltypes {
    input:
    path("data.h5ad")

    output:
    stdout

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    data = snap.read("data.h5ad", backed='r')
    labels = np.unique(data.obs['my.cell.type'])
    print('\t'.join(labels), end='')
    """
}

process pca_distance {
    input:
    tuple path("data.h5ad"), val(cell_type)

    output:
    tuple val(cell_type), path("distances.csv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import scanpy as sc
    import numpy as np
    import pandas as pd

    data = snap.read("data.h5ad", backed=None)
    data = data[data.obs['my.cell.type'] == "${cell_type}"]
    sc.pp.highly_variable_genes(data, flavor='seurat_v3', n_top_genes=3000)
    data = data[:, data.var.highly_variable]
    sc.pp.normalize_total(data, target_sum=1e4)
    sc.pp.log1p(data)
    sc.pp.scale(data, max_value=10)
    sc.pp.pca(data, n_comps=30)
    ages = sorted(np.unique(data.obs['age']))

    centroids = []
    for age in ages:
        mat = data[data.obs['age'] == age].obsm['X_pca']
        centroids.append(np.mean(mat, axis=0))
    distances = []
    for c in centroids:
        distances.append(np.linalg.norm(c - centroids[0]))
    distances = [d / max(distances) for d in distances]
    
    pd.DataFrame({
        'age': ages,
        'distance': distances
    }).to_csv('distances.csv', index=False, header=True)
    """
}

process spectral_distance {
    publishDir 'result/cell_maturation/umap', pattern: '*.pdf'

    input:
    tuple path("data.h5ad"), val(cell_type)

    output:
    tuple val(cell_type), path("distances.csv"), path("*.pdf")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import scanpy as sc
    import numpy as np
    import pandas as pd

    data = snap.read("data.h5ad", backed=None)
    data = data[data.obs['my.cell.type'] == "${cell_type}"]
    sc.pp.highly_variable_genes(data, flavor='seurat_v3', n_top_genes=3000)
    data = data[:, data.var.highly_variable]
    sc.pp.normalize_total(data, target_sum=1e4)
    sc.pp.log1p(data)
    snap.tl.spectral(data, features=None)
    snap.tl.umap(data)
    snap.pl.umap(data, color=[str(i) for i in data.obs['age']], out_file='${cell_type}.pdf')
    ages = sorted(np.unique(data.obs['age']))
    centroids = []
    for age in ages:
        mat = data[data.obs['age'] == age].obsm['X_spectral']
        centroids.append(np.mean(mat, axis=0))
    distances = []
    for c in centroids:
        distances.append(np.linalg.norm(c - centroids[0]))
    distances = [d / max(distances) for d in distances]
    pd.DataFrame({
        'age': ages,
        'distance': distances
    }).to_csv('distances.csv', index=False, header=True)
    """
}

process spectral_knn {
    input:
    tuple path("data.h5ad"), val(cell_type)

    output:
    tuple val(cell_type), path("distances.csv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import math

    data = snap.read("data.h5ad", backed=None)
    data = data[data.obs['my.cell.type'] == "${cell_type}"]
    sc.pp.highly_variable_genes(data, flavor='seurat_v3', n_top_genes=3000)
    data = data[:, data.var.highly_variable]
    sc.pp.normalize_total(data, target_sum=1e4)
    sc.pp.log1p(data)
    snap.tl.spectral(data, features=None)
    distances = snap.pp.knn(
        data,
        n_neighbors=math.ceil(math.sqrt(data.shape[0])),
        method='exact',
        inplace=False,
    )

    ages = sorted(np.unique(data.obs['age']))
    fraction = []
    matured = set(np.where(data.obs['age'] == ages[-1])[0])
    for age in ages:
        neighbors = set()
        for i in np.where(data.obs['age'] == age)[0]:
            neighbors.update(distances.getrow(i).indices)
        m = len(neighbors.intersection(matured))
        n = len(neighbors)
        fraction.append(m / n)
    fraction = [d / fraction[-1] for d in fraction]
    pd.DataFrame({
        'age': ages,
        'distance': fraction,
    }).to_csv('distances.csv', index=False, header=True)
    """
}

process combine_pca {
    publishDir "result/cell_maturation"
    input:
    tuple val(cell_type), path(distance, stageAs: "distances?.csv")

    output:
    path("*.tsv")

    """
    #!/usr/bin/env python3
    import pandas as pd

    cell_types = ${cell_type.collect {"'" + it + "'"}}
    files = ${distance.collect {"'" + it + "'"}}
    dfs = []
    for file in files:
        df = pd.read_csv(file)
        df.set_index('age', inplace=True)
        dfs.append(df.transpose())
    df = pd.concat(dfs, axis=0)
    df.insert(0, 'cell_type', cell_types)
    df.to_csv('pca_distances.tsv', index=False, header=True, sep='\t')
    """
}

process combine_spectral {
    publishDir "result/cell_maturation"
    input:
    tuple val(cell_type), path(distance, stageAs: "distances?.csv")

    output:
    path("*.tsv")

    """
    #!/usr/bin/env python3
    import pandas as pd

    cell_types = ${cell_type.collect {"'" + it + "'"}}
    files = ${distance.collect {"'" + it + "'"}}
    dfs = []
    for file in files:
        df = pd.read_csv(file)
        df.set_index('age', inplace=True)
        dfs.append(df.transpose())
    df = pd.concat(dfs, axis=0)
    df.insert(0, 'cell_type', cell_types)
    df.to_csv('spectral_distances.tsv', index=False, header=True, sep='\t')
    """
}

process combine_spectral_knn {
    publishDir "result/cell_maturation"
    input:
    tuple val(cell_type), path(distance, stageAs: "distances?.csv")

    output:
    path("*.tsv")

    """
    #!/usr/bin/env python3
    import pandas as pd

    cell_types = ${cell_type.collect {"'" + it + "'"}}
    files = ${distance.collect {"'" + it + "'"}}
    dfs = []
    for file in files:
        df = pd.read_csv(file)
        df.set_index('age', inplace=True)
        dfs.append(df.transpose())
    df = pd.concat(dfs, axis=0)
    df.insert(0, 'cell_type', cell_types)
    df.to_csv('spectral_knn.tsv', index=False, header=True, sep='\t')
    """
}

process benchmark {
    publishDir "result/cell_maturation"
    input:
    path(spectral_knn)
    path(spectral)
    path(pca)

    output:
    path("*.pdf")

    """
    #!/usr/bin/env python3
    import pandas as pd
    from plotnine import *
    from scipy.stats import spearmanr

    def cor(df):
        n, m = df.shape
        expected = list(range(m))
        return [spearmanr(df.iloc[i], expected)[0] for i in range(n)]

    spectral_knn = pd.read_csv("${spectral_knn}", sep='\t')
    spectral_knn.set_index('cell_type', inplace=True)
    spectral = pd.read_csv("${spectral}", sep='\t')
    spectral.set_index('cell_type', inplace=True)
    pca = pd.read_csv("${pca}", sep='\t')
    pca.set_index('cell_type', inplace=True)

    data = pd.DataFrame({
        'method': ['spectral_knn'] * len(spectral_knn) + ['spectral'] * len(spectral) + ['pca'] * len(pca),
        'correlation': cor(spectral_knn) + cor(spectral) + cor(pca),
    })
    p = ggplot(data, aes(x='method', y='correlation')) + geom_boxplot()
    p.save('benchmark.pdf')
    """
}

process plot_heatmap {
    publishDir "result/cell_maturation"
    input:
    path(spectral)

    output:
    path("*.pdf")

    """
    #!/usr/bin/env python3
    import pandas as pd
    from PyComplexHeatmap import *

    spectral = pd.read_csv("${spectral}", sep='\t')
    spectral.set_index('cell_type', inplace=True)

    delta = []
    for i in range(len(spectral)):
        row = spectral.iloc[i]
        delta.append([abs(row[j] - row[j + 1]) for j in range(len(row) - 1)])
    delta = pd.DataFrame(delta, index=spectral.index, columns=spectral.columns[1:])

    plt.figure(figsize=(5.5, 6.5))
    cm = ClusterMapPlotter(
        data=spectral,
        col_cluster=False,row_cluster=True,
        show_rownames=True,show_colnames=True,
        cmap='viridis',
    )
    plt.savefig("heatmap.pdf", bbox_inches='tight')

    plt.figure(figsize=(5.5, 6.5))
    cm = ClusterMapPlotter(
        data=delta,
        col_cluster=False,row_cluster=True,
        show_rownames=True,show_colnames=True,
        cmap='viridis',
    )
    plt.savefig("delta.pdf", bbox_inches='tight')
 
    """
}

workflow cell_maturation {
    take: data
    /*
    data.peak_matrix
    data.gene_matrix
    */

    main:
        cell_types = list_celltypes(data.gene_matrix) | map { it.split('\t') } | flatten

        data = data.gene_matrix | combine(cell_types)

        pca = data | pca_distance | toSortedList { a, b -> b[0] <=> a[0] } | map {it.transpose()} | combine_pca
        spectral = data | spectral_distance | map {[it[0], it[1]]} | toSortedList { a, b -> b[0] <=> a[0] } | map {it.transpose()} | combine_spectral
        spectral_knn = data | spectral_knn | toSortedList { a, b -> b[0] <=> a[0] } | map {it.transpose()} | combine_spectral_knn
        benchmark(spectral_knn, spectral, pca)

        spectral_knn | plot_heatmap
}