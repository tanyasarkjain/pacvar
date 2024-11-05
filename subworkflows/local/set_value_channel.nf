// Based on set_value_channel originally developed for isoseq
// Check input samplesheet and get read channels
//

//include { GUNZIP } from '../../modules/nf-core/gunzip/main'

workflow SET_VALUE_CHANNEL {
    take:
    infile // file: path to compressed

    main:
    Channel // Prepare value channel
          .value(file(infile))
          .set { data }


    emit:
    data // channel: [ file(infile) ]
}
