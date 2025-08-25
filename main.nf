#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// extract channels from input annotation sample sheet 
def extract_csv(csv_file) {
    // check that the sample sheet is not 1 line or less, because it'll skip all subsequent checks if so.
    file(csv_file).withReader('UTF-8') { reader ->
        def line, numberOfLinesInSampleSheet = 0;
        while ((line = reader.readLine()) != null) {
            numberOfLinesInSampleSheet++
            if (numberOfLinesInSampleSheet == 1){
                def requiredColumns = ["sample_id", "fastq_1", "fastq_2", "lane"]
                def headerColumns = line
                if (!requiredColumns.every { headerColumns.contains(it) }) {
                    log.error "Header missing or CSV file does not contain all of the required columns in the header: ${requiredColumns}"
                    System.exit(1)
                }
            }
        }
        
        if (numberOfLinesInSampleSheet < 2) {
            log.error "Provided SampleSheet has less than two lines. Provide a samplesheet with header and at least a sample."
            System.exit(1)
        }
    }

    Channel.from(csv_file)
        .splitCsv(header: true)
        .map{ row ->
            return [row.sample_id, row.fastq_1, row.fastq_2, row.lane] 
        }
}


process decontamination{
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    // container "docker://yinxiu/xenome:latest"
    container "docker://yinxiu/xengsort:latest"
    cpus 8
    memory '20 GB'
    input:
    tuple val(sample_id), path(fastq_1), path(fastq_2), val(lane)
    output:
    path("${sample_id}_l${lane}_human*graft*")
    script:
    """
    xenome classify -T ${task.cpus} -P ${params.index} --pairs --host-name ${sample_id}_mouse --graft-name ${sample_id}_l${lane}_human -i ${fastq_1} -i ${fastq_2}
    # xengsort classify --index /hpcnfs/scratch/DIMA/zhan/xengsort_index/idx --fastq ${fastq_1} --pairs ${fastq_2} --prefix ${sample_id}_l${lane}_human --mode count -T 8
    """
}

workflow {
      input=extract_csv(file(params.input))
      decontamination(input) 
}
