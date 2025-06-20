#!/usr/bin/env nextflow

/*
========================================================================================
   Variant calling Nextflow Workflow
========================================================================================
   Github   : https://github.com/scbirlab/nf-ont-call-variants
   Contact  : Eachan Johnson <eachan.johnson@crick.ac.uk>
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

/*
========================================================================================
   Help text
========================================================================================
*/
if ( params.help ) {
   println """\
         S C B I R   N A N O P O R E  V A R I A N T   C A L L I N G   P I P E L I N E
         ============================================================================
         Nextflow pipeline to call variants from Nanopore FASTQ files from bacterial clones 
         relative to a wildtype control.

         Usage:
            nextflow run scbirlab/nf-ont-call-variants --sample_sheet <csv> --inputs <dir>
            nextflow run scbirlab/nf-ont-call-variants -c <config-file>

         Required parameters:
            sample_sheet      Path to a CSV with information about the samples 
                                 and FASTQ files to be processed

         Optional parameters (with defaults):  
            inputs             Directory containing inputs. Default: "./inputs".
            trim_qual = 10     For `cutadapt`, the minimum Phred score for trimming 3' calls
            min_length = 10    For `cutadapt`, the minimum trimmed length of a read. Shorter 
                                 reads will be discarded.

         The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.
   
   """.stripIndent()
   exit 0
}

/*
========================================================================================
   Check parameters
========================================================================================
*/
if ( !params.sample_sheet ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to sample_sheet")
}

working_dir = params.outputs

processed_o = "${working_dir}/processed"
tables_o = "${working_dir}/tables"
multiqc_o = "${working_dir}/multi_qc"

log.info """\
         S C B I R   V A R I A N T   C A L L I N G   P I P E L I N E
         ===========================================================
         software
            GATK           : ${params.gatk_image}
            Clair3         : ${params.clair3_image}
            snpEff         : ${params.snpeff_url}
         inputs
            input dir.     : ${params.inputs}
            sample sheet   : ${params.sample_sheet}
         trimming 
            quality        : ${params.trim_qual}
            minimum length : ${params.min_length}
         output            
            Processed      : ${processed_o}
            Tables         : ${tables_o}
            MultiQC        : ${multiqc_o}
         """
         .stripIndent()

dirs_to_make = [processed_o, tables_o, multiqc_o]

println  """
            Making directories: 
          """.stripIndent()

dirs_to_make.each { 
   println "${it}: " 
   println file(it).mkdirs() ? "OK" : "Cannot create directory: ${it}"
}

/*
========================================================================================
   MAIN Workflow
========================================================================================
*/

workflow {

   Channel.fromPath( params.sample_sheet, 
                     checkIfExists: true )
          .splitCsv( header: true )
          .set { csv_ch }

   csv_ch
      .map { it.genome_accession }
      .unique()
      | DOWNLOAD_GENOME  // acc, fasta, gff, cds, protein

   csv_ch
      .map { tuple( 
         it.sample_id,
         it.genome_accession,
         file( 
            "${params.inputs}/${it.reads}", 
            checkIfExists: true,
         ),
      ) }
      .set { sample_ch }  // id, acc, reads

   Channel.of( [ 
      params.clair3_model, 
      params.rerio_url 
   ]) | LOAD_CLAIR3
   Channel.of( params.snpeff_url ) | LOAD_SNPEFF

   sample_ch.map { tuple( it[0], it[2] ) } | FASTQC  
   TRIM_CUTADAPT(
      sample_ch.map { tuple( it[0], it[2] ) },  // id, reads
      Channel.value( params.min_length ),
   )

   sample_ch
      .map { tuple( it[0], it[1] ) }  // id, acc
      .join( 
         TRIM_CUTADAPT.out.main,
         by: 0,
      )  // id, acc, reads
      .map { tuple( it[1], it[0], it[2] ) }  // acc, id, reads
      .set { trimmed_reads }

   DOWNLOAD_GENOME.out
      .map { tuple( it[0], it[1] ) } 
      | MINIMAP2_INDEX   // acc, minimap-idx
   trimmed_reads
      .combine( MINIMAP2_INDEX.out, by: 0 )  // acc, id, reads, minimap-idx
      .map { tuple( it[1], it[2], it[3] ) } // id, reads, minimap-idx
      | MINIMAP2_ALIGN  // id, bam, bai

   DOWNLOAD_GENOME.out
      .map { tuple( it[0], it[1] ) } 
      | BOWTIE2_INDEX   // acc, bowtie-idx-base, bowtie-idx
   trimmed_reads
      .combine( BOWTIE2_INDEX.out, by: 0 )  // acc, id, reads, bowtie-idx-base, bowtie-idx
      .map { tuple( it[1], it[2], it[3], it[4] ) }  // id, reads, bowtie-idx-base, bowtie-idx
      | BOWTIE2_ALIGN  // id, bam, bai

   DOWNLOAD_GENOME.out.map { tuple( it[0], it[1] ) } 
      | FAIDX   // acc, fai

   DOWNLOAD_GENOME.out
      .map { tuple( it[0], it[1] ) }  // acc, fasta
      .combine( FAIDX.out, by: 0 )    // acc, fasta, fai
      .set { fasta_faidx }
   
   MINIMAP2_ALIGN.out.main | MOSDEPTH 
   MINIMAP2_ALIGN.out.main | SAMTOOLS_STATS

   MINIMAP2_ALIGN.out.main     // id, bam, bai
      .combine( fasta_faidx )  // id, bam, bai, acc, fasta, fai
      .combine( LOAD_CLAIR3.out )
      | CLAIR3

   CLAIR3.out.main | NORMALIZE_VCFS

   sample_ch
      .map { tuple( it[0], it[1] ) }  // id, acc
      .join( 
         NORMALIZE_VCFS.out,
         by: 0,
      )  // id, acc, vcf, idx
      .groupTuple( by: 1 )  // [id,...], acc, [vcf, ...]
      .map { tuple( it[1], it[2] ) }  // acc, [vcf, ...]
      | MERGE_VCFS  // acc, vcf

   MERGE_VCFS.out
      .combine( DOWNLOAD_GENOME.out, by: 0 )  // acc, vcf, fasta, gff, cds, protein
      .combine( LOAD_SNPEFF.out )   // acc, vcf, fasta, gff, cds, protein, snpeff
      | SNPEFF
   SNPEFF.out.main
      | TO_TABLE

   TRIM_CUTADAPT.out.logs.concat(
         FASTQC.out.logs,
         MINIMAP2_ALIGN.out.logs,
         BOWTIE2_ALIGN.out.logs,
         SAMTOOLS_STATS.out.logs,
         SNPEFF.out.logs,
      )
      .flatten()
      .unique()
      .collect()
      | MULTIQC

}

process LOAD_SNPEFF {

   label 'some_mem'

   input:
   val snpeff_url

   output:
   path "snpEff"

   script:
   """
   curl -v -L '${snpeff_url}' > snpEff_latest_core.zip
   unzip snpEff_latest_core.zip
   """
}

process LOAD_CLAIR3 {

   tag "${clair3_model}"
   label 'some_mem'

   input:
   tuple val( clair3_model ), val( rerio_url )

   output:
   path "rerio/clair3_models/${clair3_model}"

   script:
   """
   git clone ${rerio_url}
   python rerio/download_model.py --clair3 rerio/clair3_models/${clair3_model}_model

   """
}

process DOWNLOAD_GENOME {

   tag "${accession}"
   label 'some_mem'

   input:
   val accession

   output:
   tuple val( accession ), path( "ncbi_dataset/data/*/${accession}_*_genomic.fna" ), path( "ncbi_dataset/data/*/*.gff" ), path( "ncbi_dataset/data/*/cds_from_genomic.fna" ), path( "ncbi_dataset/data/*/protein.faa" )

   script:
   """
   curl -L -A "Mozilla/5.0" --retry 3 -o genome-out "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/${accession}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"
   unzip genome-out ncbi_dataset/data/${accession}/{${accession}_*_genomic.fna,*.gff,cds_from_genomic.fna,protein.faa}
   """
}

/* 
 * Do quality control checks
 */
process FASTQC {

   label 'med_mem'

   tag "${sample_id}"

   input:
   tuple val( sample_id ), path ( reads )

   output:
   path( "*.zip" ), emit: logs

   script:
   """
   zcat ${reads} > ${sample_id}.fastq
   fastqc --noextract --memory 10000 --threads ${task.cpus} ${sample_id}.fastq
   rm ${sample_id}.fastq
   """
   stub:
   """
   zcat ${reads} | head -n1000 > ${sample_id}.fastq
   fastqc --noextract --memory 10000 --threads ${task.cpus} ${sample_id}.fastq
   rm ${sample_id}.fastq
   """

}

/*
 * Trim adapters from reads
 */
process TRIM_CUTADAPT {

   tag "${sample_id}"

   publishDir( processed_o, 
               mode: 'copy',
               pattern: "*.log" )

   input:
   tuple val( sample_id ), path( reads )
   val min_length

   output:
   tuple val( sample_id ), path( "*.trimmed.fastq.gz" ), emit: main
   path "*.log" , emit: logs

   script:
   """
   ln -s ${reads} ${sample_id}.fastq.gz
   cutadapt \
		-q ${params.trim_qual},${params.trim_qual} \
		--report=full \
      --minimum-length ${min_length} \
		-o ${sample_id}.trimmed.fastq.gz \
      --json=${sample_id}.cutadapt.json \
		${sample_id}.fastq.gz > ${sample_id}.cutadapt.log
   """
   stub:
   """
   zcat ${reads} | head -n1000 | gzip --best > ${sample_id}.fastq.gz
   cutadapt \
		-q ${params.trim_qual},${params.trim_qual} \
		--report=full \
      --minimum-length ${min_length} \
		-o ${sample_id}.trimmed.fastq.gz \
      --json=${sample_id}.cutadapt.json \
		${sample_id}.fastq.gz > ${sample_id}.cutadapt.log
   """
}

/*
 * Index the reference genome for use by Bowtie2.
 */
process BOWTIE2_INDEX {
   
   tag "${accession}"
   
   input:
   tuple val( accession ), path( genome )

   output:
   tuple val( accession ), val( "${genome.getBaseName()}" ), path( "*.bt2" )  

   script:
   """
   bowtie2-build ${genome} ${genome.getBaseName()}
   """
}

/*
 * Align reads to reference genome & create BAM file.
 */
process BOWTIE2_ALIGN {

   label 'med_mem'
   errorStrategy 'ignore'

   tag "${sample_id}"

   input:
   tuple val( sample_id ), path( reads ), val( genome_idx_base ), path( genome_idx )

   output:
   tuple val( sample_id ), path( "*.bam" ), path( "*.bai" ), emit: main
   path "*.log", emit: logs

   script:
   """
   bowtie2 \
      -x ${genome_idx_base} \
      --rdg 10,10 \
      --very-sensitive \
      --trim-to 3:2500 \
      -U ${reads} \
      2> ${sample_id}.bowtie2.log \
      | samtools sort - \
      -O bam -l 9 -o ${sample_id}.mapped.bam
   samtools index *.mapped.bam
   """
}

/*
 * Get alignment stats
 */
process MOSDEPTH {

   tag "${sample_id}"

   publishDir( path: tables_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( bamfile ), path( idx )

   output:
   path "*.txt", emit: logs
   path "*.html", emit: plots

   script:
   """
   mosdepth ${bamfile.getBaseName()} ${bamfile}
   git clone https://github.com/brentp/mosdepth.git && python mosdepth/scripts/plot-dist.py ${bamfile.getBaseName()}.mosdepth.global.dist.txt
   mv dist.html ${bamfile.getBaseName()}.mosdepth.global.dist.html
   rm -r mosdepth
   """
}

/*
 * Get alignment stats
 */
process SAMTOOLS_STATS {

   tag "${sample_id}"

   publishDir( path: tables_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( bamfile ), path( idx )

   output:
   path "*.stats", emit: logs
   path "*.png", emit: plots

   script:
   """
   samtools stats ${bamfile} > ${bamfile.getBaseName()}.stats
   plot-bamstats -p ${bamfile.getBaseName()}_plot ${bamfile.getBaseName()}.stats
   """
}

/*
 * Index the reference genome for use by Minimap2.
 */
process MINIMAP2_INDEX {

   tag "${accession}"
   
   input:
   tuple val( accession ), path( genome )

   output:
   tuple val( accession ), path( "*.mmi" )

   script:
   """
   minimap2 -x map-ont -d ${genome.getBaseName()}-ont.mmi ${genome}
   """
}

/*
 * Align reads to reference genome & create BAM file.
 */
process MINIMAP2_ALIGN {

   tag "${sample_id}" 

   publishDir( path: processed_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( reads ), path( idx )

   output:
   tuple val( sample_id ), path( "*.bam" ), path( "*.bai" ), emit: main
   path "*.log", emit: logs

   script:
   """
   minimap2 -a -c --MD -x map-ont \
      "${idx}" "${reads}" \
   > minimap.sam \
   2> "${sample_id}.minimap2.log"
   samtools sort -@ ${task.cpus} -O bam -l 9 -o "${sample_id}.mapped.bam" minimap.sam
   samtools index "${sample_id}.mapped.bam"
   rm minimap.sam
   
   """
}

/*
 * 
 */
process FAIDX {

   tag "${accession}"
   
   input:
   tuple val( accession ), path( genome )

   output:
   tuple val( accession ), path( "*.fai" ) 

   script:
   """
   samtools faidx ${genome}
   """
}


process CLAIR3 {

   tag "${sample_id}:${params.clair3_image}" 
   label 'big_cpu'

   container params.clair3_image

   input:
   tuple val( sample_id ), path( bamfile ), path( idx ), val( accession ), path( fasta ), path( fai ), path( clair3_model )

   output:
   tuple val( sample_id ), path( "merge_output.vcf.gz" ), path( fasta ), path( fai ), emit: main
   path "full_alignment.vcf.gz", emit: alignment
   path "*.log", emit: logs

   script:
   """
   set -x
   /opt/bin/run_clair3.sh \
      --bam_fn="${bamfile}" \
      --ref_fn="${fasta}" \
      --sample_name=${sample_id}\
      --threads=${task.cpus} \
      --platform="ont" \
      --model_path="${clair3_model}" \
      --output="." \
      --no_phasing_for_fa \
      --include_all_ctgs \
      --haploid_precise \
      --enable_long_indel
   ls -lah
   """
}

process NORMALIZE_VCFS {

   tag "${sample_id}" 
   label 'some_mem'

   publishDir( path: processed_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( vcf ), path( fasta ), path( fai )

   output:
   tuple val( sample_id ), path( "*.normalized.vcf.gz" )

   script:
   """
   bcftools index -t ${vcf}
   zcat ${vcf} \
      | bcftools view -i 'GT="alt"' \
      | bcftools norm -f ${fasta} -a -c e -m - \
      | bcftools norm -aD \
      | bcftools +setGT - -- -t a -n c:M \
      | bcftools sort \
      | bcftools view -i 'GT="A"' -Oz -o ${sample_id}.normalized.vcf.gz
   bcftools index --tbi ${sample_id}.normalized.vcf.gz
   """

}


process MERGE_VCFS {

   tag "${accession}"

   label 'some_mem'

   input:
   tuple val( accession ), path( vcfs )

   output:
   tuple val( accession ), path( "*-merged.vcf" )

   script:
   """
   for vcf in ${vcfs}
   do
      bcftools index -t "\$vcf"
   done
   bcftools merge --file-list <(ls -1 *.vcf.gz) -O v -o ${accession}-merged.vcf
   """
}


process SNPEFF {

   tag "${accession}:${snpeff}"

   publishDir( path: tables_o, 
               mode: 'copy' )

   input:
   tuple val( accession ), path( vcf ), path( fasta ), path( gff ), path( cds ), path( protein ), path( snpeff )

   output:
   tuple val( accession ), path( "*.{filtered,high-eff}.vcf" ), emit: main
   tuple path( "*.csv" ), path( "*.html" ), emit: logs

   script:
   """
   SNPEFF_HOME=${snpeff}
   DB_LOC=\$SNPEFF_HOME/data/${accession}
   CHR_NAME=\$(grep -v '^#' ${vcf} | cut -f1 | sort -u)

   if [ ! -e \$DB_LOC/snpEffectPredictor.bin ]
   then
      mkdir -p \$DB_LOC
      cp --remove-destination ${gff} \$DB_LOC/genes.gff
      cp --remove-destination ${fasta} \$DB_LOC/sequences.fa
      cp --remove-destination ${protein} \$DB_LOC/protein.fa
      echo '# New genome, version ${accession}' >> ${snpeff}/snpEff.config
      echo '${accession}.genome : ${accession}' >> ${snpeff}/snpEff.config
      echo '${accession}.codonTable : Bacterial_and_Plant_Plastid' >> ${snpeff}/snpEff.config
      java -jar ${snpeff}/snpEff.jar build -noCheckCds -gff3 -v ${accession}
   fi

   java -jar ${snpeff}/snpEff.jar -v -d \
      -ud 0 \
      -o gatk \
	   -csvStats ${vcf.getSimpleName()}.snpEff.csv \
      -stats ${vcf.getSimpleName()}.snpEff.html \
      ${accession} ${vcf} \
      > ${vcf.getSimpleName()}.ann.vcf

   # Remove lines where all genotypes are the same
   java -jar ${snpeff}/SnpSift.jar filter --inverse "(GEN[ALL].GT = GEN[0].GT)" \
      ${vcf.getSimpleName()}.ann.vcf \
      > ${vcf.getSimpleName()}.filtered.vcf

   # Remove silent mutants
   java -jar ${snpeff}/SnpSift.jar filter --inverse "(EFF[ALL].FUNCLASS = 'SILENT')" \
      ${vcf.getSimpleName()}.filtered.vcf \
      > ${vcf.getSimpleName()}.high-eff.vcf
   """
}

process TO_TABLE {

   tag "${accession}:${params.gatk_image}"

   publishDir( path: tables_o, 
               mode: 'copy' )

   container params.gatk_image

   input:
   tuple val( accession ), path( vcf )

   output:
   path "*.tsv"

   script:
   """
   set -x
   for f in ${vcf}
   do
      OUTFILE=\$(basename \$f .vcf).tsv
      /gatk/gatk VariantsToTable \
         -V \$f \
         -F CHROM -F POS -F REF -F ALT -F TYPE \
         -F QUAL -F FILTER -F EFF -GF GT \
         -O \$OUTFILE

      cat <(paste <(head -n1 \$OUTFILE) \
         <(printf 'Gene_Name\\tAmino_Acid_Change\\tFunctional Class')) \
         <(tail -n +2 \$OUTFILE | awk -F'|' 'BEGIN{OFS="\t"}{ print \$0,\$5,\$4,\$2 }') \
         > \$OUTFILE.1

      mv \$OUTFILE.1 \$OUTFILE
   done
   """
}

/*
 * Make log report
 */
process MULTIQC {

   publishDir( multiqc_o, 
               mode: 'copy' )

   input:
   path '*'

   output:
   tuple path( "*.html" ), path( "multiqc_data" )

   script:
   """
   multiqc .
   """
}


/*
========================================================================================
   Workflow Event Handler
========================================================================================
*/

workflow.onComplete {

   println ( workflow.success ? """
      Pipeline execution summary
      ---------------------------
      Completed at: ${workflow.complete}
      Duration    : ${workflow.duration}
      Success     : ${workflow.success}
      workDir     : ${workflow.workDir}
      exit status : ${workflow.exitStatus}
      """ : """
      Failed: ${workflow.errorReport}
      exit status : ${workflow.exitStatus}
      """
   )
}

/*
========================================================================================
   THE END
========================================================================================
*/
