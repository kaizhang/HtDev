nextflow.enable.dsl=2

include { import_data } from './pipelines/import_data';
include { atac_call_peaks } from './pipelines/atac_call_peaks';

workflow {
    data = import_data([
        "rds_file_excitatory": file(params.rna_excitatory),
        "rds_file_inhibitory": file(params.rna_inhibitory),
        "metadata_excitatory": file(params.metadata_excitatory),
        "metadata_inhibitory": file(params.metadata_inhibitory),
        "atac_fragment_files": Channel.fromPath(params.atac_files)
    ])

    atac_call_peaks([
        "merged_dataset": data.merged_atac_data
    ])
}
