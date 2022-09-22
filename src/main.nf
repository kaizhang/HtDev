nextflow.enable.dsl=2

include { import_data } from './pipelines/import_data';
include { atac_call_peaks } from './pipelines/atac_call_peaks';
include { regulatory_network } from './pipelines/regulatory_network';
include { differential_analysis } from './pipelines/differential_analysis';

workflow {
    data = import_data([
        "rds_file_excitatory": file(params.rna_excitatory),
        "rds_file_inhibitory": file(params.rna_inhibitory),
        "metadata_excitatory": file(params.metadata_excitatory),
        "metadata_inhibitory": file(params.metadata_inhibitory),
        "atac_fragment_files": Channel.fromPath(params.atac_files)
    ])

    peaks = atac_call_peaks([
        "merged_dataset": data.merged_atac_data
    ])

    regulatory_network([
        "peak_matrix": peaks.peak_matrix,
        "gene_matrix": data.merged_rna_data
    ])

    differential_analysis([
        "peak_matrix": peaks.peak_matrix,
    ])
}
