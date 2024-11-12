// Fromn the set_value_channel originally developed for isoseq
//https://github.com/nf-core/isoseq/blob/master/subworkflows/local/set_value_channel.nf

//
// Check input samplesheet and get read channels
//

include { GUNZIP } from '../../modules/nf-core/gunzip/main'

workflow SET_VALUE_CHANNEL {
    take:
    infile // file: path to compressed or not fasta/gtf

    main:
    if (infile =~ /.gz$/) {
        GUNZIP([[], file(infile)])
            .gunzip
            .map { it[1] }
            .set { data }
    }
    else {
        Channel // Prepare value channel
            .value(file(infile))
            .set { data }
    }


    emit:
    data // channel: [ file(infile) ]
}
