nextflow.enable.dsl=2

process VALIDATE_SAMPLESHEET {
    tag 'samplesheet'
    label 'small'

    input:
    path input_sheet

    output:
    path 'samplesheet.validated.csv', emit: samplesheet
    path 'samplesheet.summary.json', emit: summary

    script:
    """
    python3 ${projectDir}/bin/check_samplesheet.py \
        --input ${input_sheet} \
        --output samplesheet.validated.csv \
        --summary samplesheet.summary.json
    """
}

process FASTQC_RAW {
    tag "${meta.id}:raw"
    label 'qc'
    publishDir "${params.outdir}/qc/fastqc/raw", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*_fastqc.html'), path('*_fastqc.zip'), emit: reports

    script:
    """
    fastqc --quiet --threads ${task.cpus} ${reads}
    """
}

process TRIM_ADAPTER {
    tag "${meta.id}"
    label 'medium'
    publishDir "${params.outdir}/preprocessing/trimmed", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.cutadapt.log"), emit: logs

    script:
    """
    cutadapt \
        --cores ${task.cpus} \
        --adapter ${params.adapter} \
        --quality-cutoff ${params.trim_quality} \
        --max-n ${params.max_n_fraction} \
        --trim-n \
        --minimum-length ${Math.min(params.mirna_min_length as int, params.pirna_min_length as int)} \
        --output ${meta.id}.trimmed.fastq.gz \
        ${reads} > ${meta.id}.cutadapt.log
    """
}

process LENGTH_SPLIT {
    tag "${meta.id}"
    label 'medium'
    publishDir "${params.outdir}/preprocessing/length_split", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.miRNA.fastq.gz"), path("${meta.id}.piRNA.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.miRNA.bbduk.log"), path("${meta.id}.piRNA.bbduk.log"), emit: logs

    script:
    """
    ${params.bbduk_bin} \
        -Xmx${Math.max(1, (task.memory.toGiga() / 2) as int)}g \
        in=${reads} out=${meta.id}.miRNA.fastq.gz \
        minlength=${params.mirna_min_length} maxlength=${params.mirna_max_length} \
        minavgquality=${params.min_avg_quality} threads=${task.cpus} overwrite=t \
        2> ${meta.id}.miRNA.bbduk.log

    ${params.bbduk_bin} \
        -Xmx${Math.max(1, (task.memory.toGiga() / 2) as int)}g \
        in=${reads} out=${meta.id}.piRNA.fastq.gz \
        minlength=${params.pirna_min_length} maxlength=${params.pirna_max_length} \
        minavgquality=${params.min_avg_quality} threads=${task.cpus} overwrite=t \
        2> ${meta.id}.piRNA.bbduk.log
    """
}

process FASTQC_FILTERED {
    tag "${meta.id}:${meta.fraction}"
    label 'qc'
    publishDir "${params.outdir}/qc/fastqc/filtered", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*_fastqc.html'), path('*_fastqc.zip'), emit: reports

    script:
    """
    fastqc --quiet --threads ${task.cpus} ${reads}
    """
}

process MIRDEEP_STAR {
    tag "${meta.id}"
    label 'large'
    publishDir "${params.outdir}/mirna/mirdeep_star", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.mirdeep.counts.tsv"), emit: counts
    tuple val(meta), path("${meta.id}.mirdeep.raw.tar.gz"), emit: raw

    script:
    """
    gzip -cd ${reads} > ${meta.id}.fastq
    ln -s ${params.mirdeep_jar} MD.jar
    ln -s ${params.mirdeep_genome_dir} genome

    java -Xmx${Math.max(4, (task.memory.toGiga() - 2) as int)}g -jar MD.jar \
        -g ${params.genome} \
        -a ${params.adapter} \
        -t ${params.mirna_min_length} \
        -l ${params.mirna_max_length} \
        -s ${params.mirdeep_min_score} \
        -r ${params.mirdeep_read_cutoff} \
        -p ${params.mirdeep_parameter_p} \
        -m ${params.mirdeep_precursor_len} \
        ${meta.id}.fastq

    find . -maxdepth 2 -type f -name '*.result' -print -quit | grep -q .
    python3 ${projectDir}/bin/parse_mirdeep.py \
        --sample ${meta.id} \
        --results-dir . \
        --output ${meta.id}.mirdeep.counts.tsv

    mkdir ${meta.id}.mirdeep.raw
    find . -maxdepth 2 -type f \( -name '*.result' -o -name '*.bed' -o -name '*.html' -o -name '*.csv' \) \
        -exec cp -f {} ${meta.id}.mirdeep.raw/ \;
    tar -czf ${meta.id}.mirdeep.raw.tar.gz ${meta.id}.mirdeep.raw
    """
}

process MERGE_MIRNA_COUNTS {
    tag 'miRNA matrix'
    label 'small'
    publishDir "${params.outdir}/mirna/counts", mode: 'copy', overwrite: true

    input:
    path count_files
    path samplesheet

    output:
    path 'mirna_counts.tsv', emit: counts
    path 'mirna_candidates.tsv', emit: annotation

    script:
    """
    python3 ${projectDir}/bin/merge_counts.py \
        --inputs ${count_files.join(' ')} \
        --samplesheet ${samplesheet} \
        --counts mirna_counts.tsv \
        --annotation mirna_candidates.tsv
    """
}

process ANNOTATE_MIRNA {
    tag 'miRNA biological annotation'
    label 'medium'
    publishDir "${params.outdir}/mirna/annotation", mode: 'copy', overwrite: true

    input:
    path counts
    path candidates

    output:
    path 'mirna_counts.annotated.tsv', emit: counts
    path 'mirna_annotation.tsv', emit: annotation
    path 'novel_precursors.fasta', emit: novel_fasta

    script:
    def externalFasta = params.external_mature_fasta ? "--external-mature-fasta ${params.external_mature_fasta}" : ''
    def externalGff = params.external_mirna_gff ? "--external-gff ${params.external_mirna_gff}" : ''
    """
    python3 ${projectDir}/bin/annotate_mirna.py \
        --counts ${counts} \
        --candidates ${candidates} \
        --known-mature-fasta ${params.known_mature_fasta} \
        --known-gff ${params.known_mirna_gff} \
        ${externalFasta} ${externalGff} \
        --locus-tolerance ${params.locus_tolerance} \
        --min-samples ${params.novel_min_samples} \
        --min-total-count ${params.novel_min_total_count} \
        --min-caller-score ${params.novel_min_caller_score} \
        --output-counts mirna_counts.annotated.tsv \
        --output-annotation mirna_annotation.tsv \
        --output-fasta novel_precursors.fasta
    """
}

process RFAM_SCREEN {
    tag 'Rfam screen'
    label 'large'
    publishDir "${params.outdir}/mirna/annotation", mode: 'copy', overwrite: true

    input:
    path annotation
    path novel_fasta

    output:
    path 'mirna_annotation.rfam.tsv', emit: annotation
    path 'rfam.tblout', emit: table
    path 'rfam.log', emit: log

    script:
    """
    if [[ -s ${novel_fasta} ]]; then
        cmscan --cpu ${task.cpus} --noali --cut_ga --tblout rfam.tblout \
            ${params.rfam_cm} ${novel_fasta} > rfam.log
    else
        printf '# no novel precursor sequences were available\n' > rfam.tblout
        printf 'Rfam was skipped because the novel precursor FASTA was empty.\n' > rfam.log
    fi

    python3 ${projectDir}/bin/apply_rfam.py \
        --annotation ${annotation} \
        --tblout rfam.tblout \
        --output mirna_annotation.rfam.tsv
    """
}

process PIRNA_ALIGN_COUNT {
    tag "${meta.id}"
    label 'large'
    publishDir "${params.outdir}/pirna/alignment", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(reads), path(pirna_gff)

    output:
    tuple val(meta), path("${meta.id}.pirna.counts.tsv"), emit: counts
    tuple val(meta), path("${meta.id}.piRNA.bam"), path("${meta.id}.piRNA.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.flagstat.txt"), emit: flagstat

    script:
    def mappingMode = params.pirna_unique_only ? '-m 1' : "-k ${params.pirna_max_alignments}"
    def strandFlag = params.pirna_stranded ? '-s' : ''
    """
    gzip -cd ${reads} > ${meta.id}.piRNA.fastq
    if [[ "${pirna_gff}" == *.gz ]]; then
        gzip -cd ${pirna_gff} > piRNA.reference.gff3
    else
        cp ${pirna_gff} piRNA.reference.gff3
    fi

    ${params.bowtie_bin} -q -S --best --strata --seed ${params.bowtie_seed} \
        -p ${task.cpus} -v ${params.pirna_mismatches} ${mappingMode} --norc \
        ${params.bowtie_index} ${meta.id}.piRNA.fastq ${meta.id}.piRNA.sam

    samtools view -@ ${task.cpus} -b -F 4 ${meta.id}.piRNA.sam \
        | samtools sort -@ ${task.cpus} -o ${meta.id}.piRNA.bam
    samtools index ${meta.id}.piRNA.bam
    samtools flagstat -@ ${task.cpus} ${meta.id}.piRNA.bam > ${meta.id}.flagstat.txt

    bedtools coverage ${strandFlag} -a piRNA.reference.gff3 -b ${meta.id}.piRNA.bam -counts \
        > ${meta.id}.piRNA.coverage.gff3

    python3 ${projectDir}/bin/parse_pirna.py \
        --sample ${meta.id} \
        --input ${meta.id}.piRNA.coverage.gff3 \
        --output ${meta.id}.pirna.counts.tsv
    """
}

process MERGE_PIRNA_COUNTS {
    tag 'piRNA matrix'
    label 'small'
    publishDir "${params.outdir}/pirna/counts", mode: 'copy', overwrite: true

    input:
    path count_files
    path samplesheet

    output:
    path 'pirna_counts.tsv', emit: counts
    path 'pirna_annotation.tsv', emit: annotation

    script:
    """
    python3 ${projectDir}/bin/merge_counts.py \
        --inputs ${count_files.join(' ')} \
        --samplesheet ${samplesheet} \
        --counts pirna_counts.tsv \
        --annotation pirna_annotation.tsv
    """
}

process DESEQ2_MIRNA {
    tag 'miRNA differential expression'
    label 'medium'
    publishDir "${params.outdir}/mirna/differential_expression", mode: 'copy', overwrite: true

    input:
    path counts
    path samplesheet

    output:
    path 'miRNA_deseq2', emit: results

    script:
    """
    mkdir miRNA_deseq2
    Rscript ${projectDir}/bin/deseq2.R \
        --counts ${counts} --samplesheet ${samplesheet} --outdir miRNA_deseq2 \
        --design '${params.design}' --contrast '${params.contrast}' \
        --min-count ${params.min_count} --min-samples ${params.min_samples} --alpha ${params.alpha}
    """
}

process DESEQ2_PIRNA {
    tag 'piRNA differential expression'
    label 'medium'
    publishDir "${params.outdir}/pirna/differential_expression", mode: 'copy', overwrite: true

    input:
    path counts
    path samplesheet

    output:
    path 'piRNA_deseq2', emit: results

    script:
    """
    mkdir piRNA_deseq2
    Rscript ${projectDir}/bin/deseq2.R \
        --counts ${counts} --samplesheet ${samplesheet} --outdir piRNA_deseq2 \
        --design '${params.design}' --contrast '${params.contrast}' \
        --min-count ${params.min_count} --min-samples ${params.min_samples} --alpha ${params.alpha}
    """
}

process SOFTWARE_VERSIONS {
    tag 'versions'
    label 'small'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: true

    output:
    path 'software_versions.txt', emit: versions

    script:
    """
    {
        printf 'miRPipe-HPC\t1.0.0-alpha.1\n'
        printf 'Java\t'; java -version 2>&1 | head -n 1 || true
        printf 'FastQC\t'; fastqc --version 2>&1 | head -n 1 || true
        printf 'MultiQC\t'; multiqc --version 2>&1 | head -n 1 || true
        printf 'Cutadapt\t'; cutadapt --version 2>&1 | head -n 1 || true
        printf 'Bowtie\t'; ${params.bowtie_bin} --version 2>&1 | head -n 1 || true
        printf 'SAMtools\t'; samtools --version 2>&1 | head -n 1 || true
        printf 'BEDTools\t'; bedtools --version 2>&1 | head -n 1 || true
        printf 'R\t'; R --version 2>&1 | head -n 1 || true
        printf 'miRDeep*_jar_sha256\t'; sha256sum ${params.mirdeep_jar} 2>/dev/null | cut -d' ' -f1 || true
        printf 'container\t${params.container}\n'
    } > software_versions.txt
    """
}

process MULTIQC {
    tag 'MultiQC'
    label 'qc'
    publishDir "${params.outdir}/qc/multiqc", mode: 'copy', overwrite: true

    input:
    path qc_files

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_data', emit: data

    script:
    """
    multiqc --force --filename multiqc_report.html --outdir . .
    """
}

workflow {
    if (!params.samplesheet) {
        error 'Missing required parameter: --samplesheet'
    }
    if (!(params.run_mirna as boolean) && !(params.run_pirna as boolean)) {
        error 'At least one of --run_mirna or --run_pirna must be true.'
    }
    if (!(params.adapter ==~ /(?i)^[ACGTN]+$/)) {
        error '--adapter must contain only A, C, G, T, or N.'
    }
    if ((params.mirna_max_length as int) >= (params.pirna_min_length as int)) {
        log.warn 'miRNA and piRNA length ranges overlap. The same read can enter both branches.'
    }

    input_sheet_ch = Channel.fromPath(params.samplesheet, checkIfExists: true)
    VALIDATE_SAMPLESHEET(input_sheet_ch)

    samples_ch = VALIDATE_SAMPLESHEET.out.samplesheet
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = new LinkedHashMap(row)
            meta.id = row.sample
            meta.remove('fastq')
            tuple(meta, file(row.fastq, checkIfExists: true))
        }

    FASTQC_RAW(samples_ch)
    TRIM_ADAPTER(samples_ch)
    LENGTH_SPLIT(TRIM_ADAPTER.out.reads)

    filtered_qc_ch = LENGTH_SPLIT.out.reads.flatMap { meta, mirna_reads, pirna_reads ->
        def mirna_meta = new LinkedHashMap(meta)
        def pirna_meta = new LinkedHashMap(meta)
        mirna_meta.fraction = 'miRNA'
        pirna_meta.fraction = 'piRNA'
        [tuple(mirna_meta, mirna_reads), tuple(pirna_meta, pirna_reads)]
    }
    FASTQC_FILTERED(filtered_qc_ch)

    qc_files_ch = FASTQC_RAW.out.reports.flatMap { meta, html, zip -> [html, zip] }
        .mix(TRIM_ADAPTER.out.logs.map { meta, log_file -> log_file })
        .mix(LENGTH_SPLIT.out.logs.flatMap { meta, mirna_log, pirna_log -> [mirna_log, pirna_log] })
        .mix(FASTQC_FILTERED.out.reports.flatMap { meta, html, zip -> [html, zip] })

    if (params.run_mirna as boolean) {
        mirna_reads_ch = LENGTH_SPLIT.out.reads.map { meta, mirna_reads, pirna_reads -> tuple(meta, mirna_reads) }
        MIRDEEP_STAR(mirna_reads_ch)
        MERGE_MIRNA_COUNTS(MIRDEEP_STAR.out.counts.map { meta, counts -> counts }.collect(), VALIDATE_SAMPLESHEET.out.samplesheet)
        ANNOTATE_MIRNA(MERGE_MIRNA_COUNTS.out.counts, MERGE_MIRNA_COUNTS.out.annotation)

        if (params.run_rfam as boolean) {
            RFAM_SCREEN(ANNOTATE_MIRNA.out.annotation, ANNOTATE_MIRNA.out.novel_fasta)
            qc_files_ch = qc_files_ch.mix(RFAM_SCREEN.out.log)
        }

        if (params.run_deseq2 as boolean) {
            DESEQ2_MIRNA(ANNOTATE_MIRNA.out.counts, VALIDATE_SAMPLESHEET.out.samplesheet)
        }
    }

    if (params.run_pirna as boolean) {
        pirna_gff_path = params.pirna_gff ?: "${projectDir}/../refs/pirnadb/pirnadb.${params.genome}.gff3.gz"
        pirna_gff = file(pirna_gff_path, checkIfExists: true)
        pirna_reads_ch = LENGTH_SPLIT.out.reads.map { meta, mirna_reads, pirna_reads -> tuple(meta, pirna_reads, pirna_gff) }
        PIRNA_ALIGN_COUNT(pirna_reads_ch)
        MERGE_PIRNA_COUNTS(PIRNA_ALIGN_COUNT.out.counts.map { meta, counts -> counts }.collect(), VALIDATE_SAMPLESHEET.out.samplesheet)
        qc_files_ch = qc_files_ch.mix(PIRNA_ALIGN_COUNT.out.flagstat.map { meta, flagstat -> flagstat })

        if (params.run_deseq2 as boolean) {
            DESEQ2_PIRNA(MERGE_PIRNA_COUNTS.out.counts, VALIDATE_SAMPLESHEET.out.samplesheet)
        }
    }

    SOFTWARE_VERSIONS()
    MULTIQC(qc_files_ch.mix(SOFTWARE_VERSIONS.out.versions).collect())
}
