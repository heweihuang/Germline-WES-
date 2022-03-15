#!/bin/env nextflow

/*
 * Releases
 */

// user: heweihuang
// version: v1.0.0  whhe(2022-01-21)

/*
 * print usage
 */
params.h = false
params.help = false
if(params.help || params.h){
log.info ''
log.info 'ST(somatic tumor)_DNA.nf ~ somatic tumor DNA Pipeline v1.0.0'
log.info '=================================================================================================================================='
log.info 'call the var of germline pipeline for rawdata fastq.gz'
log.info 'Usage:'
log.info 'germline_1.0.0.nf --i /rawdata input path [*fq.gz]/ --t cpu_number --o /outut path/ '
log.info ''
log.info 'Options:'
log.info '      --help                      Show this message and exit.'
log.info '      --i                 <str>   input (fastq.gz/fq.gz) path                              [./]'
log.info '      --id                <str>   the sample name                                          [out]'
log.info '      --o                 <str>   the output dir                                           [./]'
log.info '      --normal            <str>   the Normal sample bam                                    []'
log.info '      --panel_type        <str>   panel type'
log.info '      --groups_type       <str>   groups type'
log.info '  QC_options:'
log.info '      --qc                <str>   Whether the rawdata is qc filtered                       [true]'
log.info '      --qc_cpu            <int>   cpu number for each QC job                               [4]'
log.info '      --adapter           <str>   adapter sequences file                                   [/fuer2/01.software/2.fastp/adapters/TruSeq3-PE-2.fa]'
log.info '      --clean_rate        <int>   threshold of clean data rate                             [80]'
log.info '  Map_options:'
log.info '      --m                 <str>   Whether the reads is mapping the reference database      [true]'
log.info '      --map_cpu           <int>   cpu number for each rm rRNA job                          [8]'
log.info '      --mapping_rate      <int>   threshold of rRNA mapping rate                           [80]'
log.info '      --target_bed        <str>   Capture range of the panel                               [/fuer2/01.Pipeline/03.pipeline_test/01.ST_DNA/db/bed/regioncount.bed]'
log.info '      --expan_bed         <str>   panel bed expan range                                    [/fuer2/01.Pipeline/03.pipeline_test/01.ST_DNA/db/bed/regioncount.expan.bed]'
log.info '  SNV_options:'
log.info '      --c                 <str>   Whether the bam is call snp and indel variation          [true]'
log.info '      --snv_cpu           <int>   cpu number for each call snp and indel job               [1]'
log.info '      --scatter_count     <int>   How many parts of the bed divided for parallel analysis  [default/10]'
log.info '      --list              <str>   split list of mutect call '
log.info '      --d                 <str>   the bed is divided into multiple parts to analysis       [true]'
log.info '  SV_options:'
log.info '      --s                 <str>   Whether the bam is call sv variation                     [true]'
log.info '      --sv_cpu            <int>   cpu number for each call sv job                          [8]'
log.info '  CNV_options:'
log.info '      --c                 <str>   Whether the bam is call cnv                              [true]'
log.info '      --cnv_cpu           <int>   cpu number for each call cnv job                         [8]'
log.info '  Fusion_options:'
log.info '      --f                 <str>   Whether the bam is call fusion                           [true]'
log.info '      --cnv_cpu           <int>   cpu number for each call cnv job                         [8]'
log.info '  UMI_options:'
log.info '      --u                 <str>   Whether the bam is call UMI                              [true]'
log.info '      --umi_cpu           <int>   cpu number for each call umi job                         [8]'
log.info '  Gemr_options:'
log.info '      --g                 <str>   Whether the bam is call germline                         [true]'
log.info '      --gemr_cpu          <int>   cpu number for each call gemr job                        [10]'
log.info '   Anno_options:'
log.info '      --a                 <str>   Whether the vcf is annotation DB                         [true]'               
log.info '      --anno_cpu          <int>   cpu number for each DB anno job                          [4]'
log.info '   Result:'
log.info '   --dm                   <str>   data merge to report                                     [true]'
exit 1
}

outdirAbs              = ""

if( params.d[0] == "." ){
    outdirAbs = "${launchDir}/${params.d}"
}else{
    outdirAbs = "${params.d}"
}

if(params.project_type == "1"){
    Tbed    = "${params.target_bed1}"
    Tbed500 = "${params.expand_bed1}"
       
  }
if(params.project_type == "2"){
    Tbed    = "${params.target_bed2}"
    Tbed500 = "${params.expand_bed2}" 
  }

Channel
    .fromFilePairs("$params.i/*_{1,2}.{fq,fastq}*",size:2)
    .ifEmpty { error "Cannot find fq/fq.gz/fastq/fastq.gz files: $params.i" }
    .set{raw_reads}

Channel
    .fromPath("${params.list}/*.list")
    .flatten()
    .set{intervals}

Channel
    .fromPath("${params.list}/*.list")
    .flatten()
    .set{intervals1}

Channel
    .fromPath("${params.list}/*.list")
    .flatten()
    .set{intervals2}

Channel
    .fromPath("${params.list}/*.list")
    .flatten()
    .set{intervals3}

Channel
    .fromPath("${params.list}/*.list")
    .flatten()
    .set{intervals4}

process QC{
    tag "${libName}"
    cpus "${params.qc_cpu}"
   
    input:
    set val(libNameRaw), file(r_reads) from raw_reads
    output:
    set val(libName), file("*Clean.fastq.gz") into clean_fq1, clean_fq2, clean_fq3, clean_fq4, clean_fq5, clean_fq6, clean_fq7, clean_fq8, clean_fq9
    file("*fastp.html") into fastp_log
    file("*rawDataStat.xls") into fastp_raw_stat
    set val(libName), file("*filter_ratio.xls") into fastp_clean_stat, fastp_clean_stat2
    file("*md5.txt")
    script:
    if(params.groups_type != "6"){ 
       libName = libNameRaw.split(/[-_]/)[0] + "-" + libNameRaw.split(/[-_]/)[1]
    }else{
       libName = libNameRaw.split(/[-_]/)[2] + "-" + libNameRaw.split(/[-_]/)[3] + "-" + libNameRaw.split(/[-_]/)[6]

    }

    """
    ${params.fastp} --thread ${params.qc_cpu} -i ${r_reads[0]} -I ${r_reads[1]} -o ${libName}.read1_Clean.fastq.gz -O ${libName}.read2_Clean.fastq.gz --adapter_fasta=${params.adapter} --html=${libName}.fastp.html 1>fastp.log 2>fastp.err
    python3 ${params.fastp_stat} fastp.err ./ ${libName}
    md5sum *Clean.fastq.gz > ${libName}.md5.txt
    """ 
}   

process check_QC {
    tag "${libName}"
    cpus "1"
    errorStrategy='finish'

    input:
    set val(libName), file(qc_s) from fastp_clean_stat
    """
    perl ${params.check_cleanR} -i ${qc_s} -threshold ${params.clean_rate}
    """
}

if( params.UMI == "1"){
process extra_UMI {
    tag "${libName}"
    cpus "${params.UMI_cpu}"

    input:
    set val(libName), file(c_reads) from clean_fq8

    output:
    set val(libName), file("*.markdup.sort.bam"), file("*.markdup.sort.bam.bai") into mapping_bam1, mapping_bam2, mapping_bam3, mapping_bam4, mapping_bam5, mapping_bam6, mapping_bam7,mapping_bam8, mapping_bam9, mapping_bam10, mapping_bam11,mapping_bam12, mapping_bam13, mapping_bam14, mapping_bam15, mapping_bam16, mapping_bam17, mapping_bam18, mapping_bam19, mapping_bam20, mapping_bam21, mapping_bam22
    set val(libName), file("*.markdup.sort.bam") into mp_bam, mp_bam1, mp_bam2, mp_bam3, mp_bam4, mp_bam5, mp_bam6, mp_bam7
    """
    java -Xmx8G -jar ${params.picard} FastqToSam F1=${c_reads[0]} F2=${c_reads[1]} O=${libName}.unmapped.bam SAMPLE_NAME=${libName}
    java -jar ${params.fgbio} ExtractUmisFromBam --input=${libName}.unmapped.bam --output=${libName}.unmapped.withUMI.bam --read-structure=10M140T 10M140T --molecular-index-tags=ZA ZB --single-tag=RX
    java -Xmx8G -jar ${params.picard} SamToFastq I=${libName}.unmapped.withUMI.bam F=${libName}.SamToFastq INTERLEAVE=true
    ${params.bwa} mem -t ${params.UMI_cpu} -p ${params.fa} ${libName}.SamToFastq > ${libName}.mapped.bam
    java -Xmx8G -jar ${params.picard} MergeBamAlignment UNMAPPED=${libName}.unmapped.withUMI.bam ALIGNED=${libName}.mapped.bam O=${libName}.mapped.withUMI.bam R=${params.fa} SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true
    java -Xmx8G -jar ${params.fgbio} GroupReadsByUmi --input=${libName}.mapped.withUMI.bam --output=${libName}.GroupedReads.bam --strategy=paired --edits=1 --min-map-q=20
    java -Xmx8G -jar ${params.fgbio} CallDuplexConsensusReads --input=${libName}.GroupedReads.bam --output=${libName}.consensus.unmapped.bam --error-rate-pre-umi=45 --error-rate-post-umi=30 --min-input-base-quality=30
    java -Xmx8G -jar ${params.fgbio} FilterConsensusReads --input=${libName}.consensus.unmapped.bam --output=${libName}.consensus.filtered.unmapped.bam --ref=${params.fa} --min-reads=2 1 1 --max-read-error-rate=0.05 --max-base-error-rate=0.1 --min-base-quality=50 --max-no-call-fraction=0.05
    java -jar ${params.picard} SamToFastq I=${libName}.consensus.filtered.unmapped.bam F=${libName}.consensus.filtered.unmapped.Fastq INTERLEAVE=true
    ${params.bwa} mem -t ${params.UMI_cpu} -p ${params.fa} ${libName}.consensus.filtered.unmapped.Fastq  >  ${libName}.consensus.filtered.mapped.bam
    java -Xmx8G -jar ${params.picard} SortSam INPUT=${libName}.consensus.filtered.mapped.bam   OUTPUT=${libName}.consensus.filtered.mapped.sorted.bam  SORT_ORDER=queryname
    java -Xmx8G -jar ${params.picard} SortSam INPUT=${libName}.consensus.filtered.unmapped.bam   OUTPUT=${libName}.consensus.filtered.unmapped.sorted.bam  SORT_ORDER=queryname
    java -Xmx8G -jar ${params.picard} MergeBamAlignment UNMAPPED=${libName}.consensus.filtered.unmapped.sorted.bam  ALIGNED=${libName}.consensus.filtered.mapped.sorted.bam O=${libName}.consensus.filtered.mapped.withUMI.bam R=${params.fa} SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true
    java -Xmx8G -jar ${params.fgbio} ClipBam --input=${libName}.consensus.filtered.mapped.withUMI.bam --output=${libName}.consensus.filtered.mapped.withUMI.clipped.bam --ref=${params.fa} --soft-clip=false --clip-overlapping-reads=true
    ln -s ${libName}.consensus.filtered.mapped.withUMI.clipped.bam ${libName}.markdup.sort.bam
    ${params.samtools} index ${libName}.markdup.sort.bam
    ln -s ${libName}.consensus.filtered.mapped.withUMI.clipped.bam ${libName}.bqsr.bam
    ln -s ${libName}.markdup.sort.bam.bai ${libName}.bqsr.bam.bai
    """
 }
}else{
process map_bwa {
    tag "${libName}"
    cpus "${params.mapping_cpu}"

    input:
    set val(libName), file(c_reads1) from clean_fq1
    
    output:
    set val(libName), file("*.markdup.sort.bam"), file("*.markdup.sort.bam.bai") into mapping_bam1, mapping_bam2, mapping_bam3, mapping_bam4, mapping_bam5, mapping_bam6, mapping_bam7,mapping_bam8, mapping_bam9, mapping_bam10, mapping_bam11,mapping_bam12, mapping_bam13, mapping_bam14, mapping_bam15, mapping_bam16, mapping_bam17, mapping_bam18, mapping_bam19, mapping_bam20, mapping_bam21, mapping_bam22, mapping_bam23
    set val(libName), file("*.markdup.sort.bam") into mp_bam, mp_bam1, mp_bam2, mp_bam3, mp_bam4, mp_bam5, mp_bam6, mp_bam7
    set val(libName), file("*.chromosome.stat") into chr_stat, chr_stat1
    """
    ${params.bwa} mem -t ${params.mapping_cpu} -R "@RG\\tID:${libName}\\tSM:${libName}\\tPL:illumina\\tCN:BDLS\\tLB:BDLS" ${params.fa} ${c_reads1[0]} ${c_reads1[1]} | ${params.samtools} view -1 -bS - -o ${libName}.bam 
    ${params.gatk} SortSam -I ${libName}.bam -O ${libName}.sort.bam -SO coordinate
    ${params.gatk} MarkDuplicates -I ${libName}.sort.bam -O ${libName}.markdup.sort.bam -M markdup.sort.bam.txt 
    ${params.samtools} index ${libName}.markdup.sort.bam
    ${params.samtools} view ${libName}.markdup.sort.bam|perl -lane '{\$h{\$F[2]}++}END{print "\$_\\t\$h{\$_}" foreach sort keys %h}' > ${libName}.chromosome.stat
    """
  }
}


process IGV_url{
    tag "${libName}"
    cpus "2"

    input:
    set val(libName), file(bam), file(bai) from mapping_bam23
    """
    /fuer2/01.Pipeline/03.pipeline_test/whhe/01.ST_DNA/scr/rsync_auto.ex ${libName}.markdup.sort.bam* /cinomedical/nfsdata/02.software/jbrowse/data/hg19/bam/ 
    """  
}


process dep_cov{
    tag "${libName}"
    cpus "2"
    publishDir "${params.o}/01.QC/${libName}/", mode: 'copy', overwrite: true    
    input:
    set val(libName), file(f) from mp_bam3.concat(fastp_clean_stat2).groupTuple()
    output: 
    set val(libName), file("*depth.tsv.gz") into depth_file
    set val(libName), file("*.qcStat.xls")  into qc_result
    """
    python3 ${params.count_dep} -o dep_cov -b ${Tbed} -bam ${f[0]} -t 1,10,20,50,100,200,500,1000 
   # perl ${params.map_stat} ./dep_cov ${libName} ${libName}.quality.xls   
    ln -s ./dep_cov/*depth.tsv.gz ${libName}.depth.tsv.gz
    perl ${params.dep_cov_stat} -n ${libName} -p ${libName}.quality -gb ${params.hg19_gene_bed} -ucb dep_cov/uncover.bed dep_cov/coverage.report dep_cov/cov_summary.txt
    python3 ${params.qc_stat} --qc ${f[1]} --map ${libName}.quality.xls --o ${libName}.qcStat.tmp.xls
    awk -F '\\t' 'NR==2||NR==3{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS \$i:\$i}END{for(i=0;i++<NF;)print a[i]}' dep_cov/cov_summary.txt|paste ${libName}.qcStat.tmp.xls - > ${libName}.qcStat.xls
    """   
}

process BQSR{
    tag "${libName}"
    cpus "1"
  
    input:
    set val(libName), file(bam), file(bai) from mapping_bam2
    output:
    set val(libName), file("*.bqsr.bam"), file("*.bqsr.bam.bai") into bqsr_bam1, bqsr_bam2, bqsr_bam3, bqsr_bam4, bqsr_bam5, bqsr_bam6
    """
    ${params.gatk} --java-options -Xmx20000m BaseRecalibrator -I $bam -R ${params.fa} --known-sites ${params.dbsnp138} --known-sites ${params.indel_phase1} --known-sites ${params.mills} -O ${libName}_recal.table
    ${params.gatk} ApplyBQSR -R ${params.fa} -I $bam --bqsr-recal-file ${libName}_recal.table -O ${libName}.bqsr.bam
    ${params.sambamba} index ${libName}.bqsr.bam 
    """
}

if(params.scatter_count != "default" && params.scatter_count != "0"){
process SplitIntervals{
    tag "${libName}"
    cpus "1" 
   
    input: 
    file(${Tbed}) 
    output:
    file("*.interval_list") into interval_List    
    """
    ${params.gatk} --java-options -Xmx3000m SplitIntervals -R ${params.fa} -L ${Tbed} -scatter ${params.scatter_count} -O ${params.o}/03.SNV/${libName}/interval-files 
    """       
    }
}


if(params.scatter_count != "default" && params.scatter_count != "0"){
process M2_scatter{
    tag "${libName}_${bed[0]}_hc"
    cpus "1"
 
    input:
          set file(bed), val(libName), file(bam) from intervals_list.combine(bqsr_bam1)     
    output:
          set val(libName), file("*.g.vcf"), file("*vcf.idx") into HC_result
    """
    NAME=`ls ${bed[0]} | awk -F '.' '{print \$1i"."\$2}'
    ${params.gatk} HaplotypeCaller -I ${bam} -R ${params.fa} -O ${libName}.germline.\${NAME}.vcf -L ${bed} -ERC GVCF --native-pair-hmm-threads ${params.HC_cpu}
    """
  }
}else{


process HC_bed{
    tag "${libName}_${bed[0]}_hc}"
    cpus "${params.HC_cpu}"

    input:
        set file(bed), val(libName), file(bam),file(bai) from intervals1.combine(bqsr_bam2)

    output:
    set val(libName), file("*.vcf"), file("*vcf.idx") into HC_result1, HC_result2
    """
    NAME=`ls ${bed} | awk -F '.' '{print \$1i"."\$2}'`
    ${params.gatk} HaplotypeCaller -I ${bam} -R ${params.fa} -O ${libName}.germline.\${NAME}.vcf -L ${bed} -ERC GVCF --native-pair-hmm-threads ${params.HC_cpu}
    """
   }
}

process HC_merge{
    tag "${libName}"
    cpus "1"
    
    input:
    set val(libName), file(vcf), file(idx) from HC_result1.groupTuple(by: 0)
    output:
    set val(libName), file("*merge.vcf") into vcf_merge, vcf_merge1
    """
    VCFs=`ls *vcf |perl -ne 'chomp;print "-I \$_ "'`
    ${params.gatk} --java-options -Xmx4000m MergeVcfs \${VCFs} -O ${libName}.merge.vcf
    """
}


if(params.groups_type == "1"){
process merge{
    tag "${libName}"
    cpus "1"

    input:
    set val(libName), file(vcf), file(idx) from vcf_merge.collect()
    output:
    set val(libName), file("*GenotypeGVCFs.snp.vcf") into snp_vcf
    set val(libName), file("*GenotypeGVCFs.indel.vcf") into indel_vcf
    """
    VCFs=`ls *vcf |perl -ne 'chomp;print "-V \$_ "' `
    ${params.gatk} CombineGVCFs \${VCFs} -R ${params.fa} -O ${libName}.merge.g.vcf
    ${params.gatk} GenotypeGVCFs -R ${params.fa} -V ${libName}.merge.g.vcf -O ${libName}.GenotypeGVCFs.vcf -L ${Tbed}
    ${params.gatk} SelectVariants -select-type SNP -V ${libName}.GenotypeGVCFs.vcf -O ${libName}.GenotypeGVCFs.snp.vcf
    ${params.gatk} SelectVariants -select-type INDEL -V ${libName}.GenotypeGVCFs.vcf -O ${libName}.GenotypeGVCFs.indel.vcf
    """
    }
}else{
process trio_merge{
    tag "${params.kid}"
    cpus "1"

    input:
    set val(libName1), file(vcf1), val(libName2), file(vcf2), val(libName3), file(vcf3) from vcf_merge.collect()
    output:
    set val(params.kid), file("*GenotypeGVCFs.snp.vcf") into snp_vcf
    set val(params.kid), file("*GenotypeGVCFs.indel.vcf") into indel_vcf
    """
    VCFs=`ls *vcf |perl -ne 'chomp;print "-V \$_ "' `
    ${params.gatk} CombineGVCFs \${VCFs} -R ${params.fa} -O ${params.kid}.merge.g.vcf
    ${params.gatk} GenotypeGVCFs -R ${params.fa} -V ${params.kid}.merge.g.vcf -O ${params.kid}.GenotypeGVCFs.vcf -L ${Tbed}
    ${params.gatk} SelectVariants -select-type SNP -V ${params.kid}.GenotypeGVCFs.vcf -O ${params.kid}.GenotypeGVCFs.snp.vcf
    ${params.gatk} SelectVariants -select-type INDEL -V ${params.kid}.GenotypeGVCFs.vcf -O ${params.kid}.GenotypeGVCFs.indel.vcf
    """
    }
}

process snp_filter{
    tag "${libName}"
    cpus "1"
    input:
    set val(libName) ,file(snp) from snp_vcf
    output:
    set val(libName), file("*GenotypeGVCFs.VQSR.Merge.snp.filter.vcf") into f_snp_vcf
    """
    ${params.gatk} VariantRecalibrator -R ${params.fa} -V ${snp} --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.hapmap} --resource:omni,known=false,training=true,truth=false,prior=12.0 ${params.G1000_omni} -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${params.G1000_phase1} -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.dbsnp_138} -mode SNP -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 --rscript-file ${libName}.VQSR.snp.plots.R --tranches-file ${libName}.snp.tranches -O ${libName}.VQSR.snp.recal --max-gaussians 4
    ${params.gatk} ApplyVQSR -R ${params.fa} -V ${snp} --recal-file ${libName}.VQSR.snp.recal --mode SNP -ts-filter-level 99.0 --tranches-file ${libName}.snp.tranches -O ${libName}.GenotypeGVCFs.VQSR.snp.vcf
    ${params.gatk} VariantFiltration --variant ${libName}.GenotypeGVCFs.VQSR.snp.vcf -O ${libName}.GenotypeGVCFs.VQSR.Merge.snp.filter.vcf --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filter-name "Standard" --filter-expression "DP <= 5" --filter-name "DP_5" --filter-expression  "DP < 20" --filter-name "DP_20"
    """      
}

process indel_filter{
    tag "${libName}"
    cpus "1"
    input:
    set val(libName), file(indel) from indel_vcf
    output:
    set val(libName), file("*GenotypeGVCFs.VQSR.Merge.indel.filter.vcf") into f_indel_vcf
    """
    ${params.gatk} VariantRecalibrator -R ${params.fa} -V ${indel} -resource:mills,known=false,training=true,truth=true,prior=12.0 ${params.mills} -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.dbsnp_138} -mode INDEL -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum --rscript-file ${libName}.VQSR.indel.plots.R --tranches-file ${libName}.VQSR.indel.tranches -O ${libName}.VQSR.indel.recal --max-gaussians 2
    ${params.gatk} ApplyVQSR -R ${params.fa} -V ${indel} --recal-file ${libName}.VQSR.indel.recal --mode INDEL --truth-sensitivity-filter-level 99.0 --tranches-file ${libName}.VQSR.indel.tranches -O ${libName}.GenotypeGVCFs.VQSR.indel.vcf
    ${params.gatk} VariantFiltration --variant ${libName}.GenotypeGVCFs.VQSR.indel.vcf -O ${libName}.GenotypeGVCFs.VQSR.Merge.indel.filter.vcf --filter-expression "QD < 2.0 || ReadPosRankSum < -20.0 || FS > 200.0 || SOR > 10.0" --filter-name "Standard" --filter-expression "DP <= 5" --filter-name "DP_5" --filter-expression  "DP < 20" --filter-name "DP_20"
    """
}

process merge_vcf{
    tag "${libName}"
    cpus "1"
    input:
    set val(libName), file(args) from f_snp_vcf.mix(f_indel_vcf).groupTuple(by: 0)
    output:
    set val(libName), file("*.merge.filter.vcf") into filter_vcf, filter_vcf1
    """
    ${params.gatk} MergeVcfs -I ${libName}.GenotypeGVCFs.VQSR.Merge.snp.filter.vcf -I ${libName}.GenotypeGVCFs.VQSR.Merge.indel.filter.vcf -O ${libName}.merge.filter.vcf
    """
}

if(params.groups_type != "1"){
    process trio {
    tag "${params.kid}"
    cpus "1"
    //errorStrategy='ignore'
    input:
    set val(params.son), file(vcf) from filter_vcf
    output:
    set val(params.son), file("*pass.trio.vcf") into trio_vcf
    """
    echo -e "FAM1\\t${params.father}\\t0\\t0\\t1\\t${params.f_ill}\\nFAM1\\t${params.mother}\\t0\\t0\\t2\\t${params.m_ill}\\nFAM1\\t${params.kid}\\t${params.father}\\t${params.mother}\\t${params.kid_sex}\\t${params.kid_ill}" > fam.ped
    awk -F '\\t' -vs=${params.kid} '{if(\$0~/proband/){print "    proband: "s}else{print \$0}}' ${params.exomiser_cfg} > analysis-exome.yml
    ln -s ${params.application}
    mkdir ./results
    CPATH=${params.exomiser_cfg}
    perl ${params.vcf_rmchr} ${vcf} > in.vcf
    ssh mgr "cd `pwd`; /fuer2/03.Soft/01.Soft_project/java/jdk-17.0.2/bin/java -jar ${params.exomiser} --analysis analysis-exome.yml"
    grep "#" results/sample-trio-FULL_AD.vcf > results/head.txt
    awk  '\$0!~/#/' results/*.vcf|sort -V|cat results/head.txt - > results/sample.trio.vcf
    awk '\$7~/PASS/||\$1~/#/' results/sample.trio.vcf  > results/sample.pass.trio.vcf
    ln -s results/sample.pass.trio.vcf ${params.kid}.pass.trio.vcf
    """
  }
process trio_anno {
    tag "${params.kid}"
    cpus "1"

    input:
    set val(libName), file(fvcf) from trio_vcf
    output:
    set val(libName), file("*anno.hg19_multianno.final.clinvar_hgmd.txt") into anno_vcf1

    """
    ${params.convert2annovar} -format vcf4old --includeinfo ${fvcf} > ${libName}.variant
    ${params.table_annovar} ${libName}.variant -thread 16 ${params.anno_db} --buildver hg19 -remove -protocol refGene,dbnsfp31a_interpro,avsnp150,1000g2015aug_all,exac03,gnomad_genome,intervar_20180118,clinvar_20200316,HGMD_generic_short,cosmic91_revised,dbscsnv11,dbnsfp35c,genomicSuperDups,gwascatalog,ljb26_all,omim201806 -operation gx,f,f,f,f,f,f,f,f,f,f,f,r,r,f,r -arg '-hgvs',,,,,,,,,,,,,,, -otherinfo -nastring . -xreffile ${params.anno_db}/gene_fullxref.txt -out ${libName}.anno
    awk -F '\\t' '\$177!~/Standard/ && \$177!~/DP_5/' ${libName}.anno.hg19_multianno.txt|cut -f 1-10,17,26-29,32,37,41,45,74-141,167-169,180 > ${libName}.anno.hg19_multianno.filter1.txt
    awk -F '\\t' '(NR==1) || ((\$6=="exonic" || \$6=="splicing") && (\$9!~/^synonymou/ && \$9!~/^unknown/) && ((\$16!="." && \$16<0.01)||(\$16=="." && \$15!="." && \$15<0.01 )||(\$16=="." && \$15=="." && \$14!="." && \$14<0.01 )||(\$16=="." && \$15=="." && \$14=="." && \$17<0.01)))' ${libName}.anno.hg19_multianno.filter1.txt|grep -vEw "strand_bias|normal_artifact|strand_bias|weak_evidence|base_qual|multiallelic|base_qual|fragment|position|map_qual|slippage|synonymous|panel_of_normals" > ${libName}.anno.hg19_multianno.filter2.txt
    cat ${libName}.anno.hg19_multianno.filter2.txt | sed 's/Otherinfo/GT\tphredscore\tDP/'|cut -f1-90 > ${libName}.anno.hg19_multianno.filter3.txt
    cut -f 91 ${libName}.anno.hg19_multianno.filter2.txt | cut -d: -f 2 | sed '1d' | awk -F ',' 'BEGIN{print "AF"}(\$1+\$2)!=0{print \$2/(\$1+\$2)}(\$1+\$2)==0{print "both is 0"}'| paste -d "\t"  ${libName}.anno.hg19_multianno.filter3.txt - > ${libName}.anno.hg19_multianno.final.txt
    perl ${params.filter_ch} ${libName}.anno.hg19_multianno.final.txt > ${libName}.anno.hg19_multianno.final.clinvar_hgmd.txt
    """
 }
}else{
process anno {
    tag "${libName}"
    cpus "1"

    input:
    set val(libName), file(fvcf) from filter_vcf1
    output:
    set val(libName), file("*anno.hg19_multianno.final.clinvar_hgmd.txt") into anno_vcf1

    """
    ${params.convert2annovar} -format vcf4old --includeinfo ${fvcf} > ${libName}.variant
    ${params.table_annovar} ${libName}.variant -thread 16 ${params.anno_db} --buildver hg19 -remove -protocol refGene,dbnsfp31a_interpro,avsnp150,1000g2015aug_all,exac03,gnomad_genome,intervar_20180118,clinvar_20200316,HGMD_generic_short,cosmic91_revised,dbscsnv11,dbnsfp35c,genomicSuperDups,gwascatalog,ljb26_all,omim201806 -operation gx,f,f,f,f,f,f,f,f,f,f,f,r,r,f,r -arg '-hgvs',,,,,,,,,,,,,,, -otherinfo -nastring . -xreffile ${params.anno_db}/gene_fullxref.txt -out ${libName}.anno
    awk -F '\\t' '\$177!~/Standard/ && \$177!~/DP_5/' ${libName}.anno.hg19_multianno.txt|cut -f 1-10,17,26-29,32,37,41,45,74-141,167-169,180 > ${libName}.anno.hg19_multianno.filter1.txt
    awk -F '\\t' '(NR==1) || ((\$6=="exonic" || \$6=="splicing") && (\$9!~/^synonymou/ && \$9!~/^unknown/) && ((\$16!="." && \$16<0.01)||(\$16=="." && \$15!="." && \$15<0.01 )||(\$16=="." && \$15=="." && \$14!="." && \$14<0.01 )||(\$16=="." && \$15=="." && \$14=="." && \$17<0.01)))' ${libName}.anno.hg19_multianno.filter1.txt|grep -vEw "strand_bias|normal_artifact|strand_bias|weak_evidence|base_qual|multiallelic|base_qual|fragment|position|map_qual|slippage|synonymous|panel_of_normals" > ${libName}.anno.hg19_multianno.filter2.txt
    cat ${libName}.anno.hg19_multianno.filter2.txt | sed 's/Otherinfo/GT\tphredscore\tDP/'|cut -f1-90 > ${libName}.anno.hg19_multianno.filter3.txt
    cut -f 91 ${libName}.anno.hg19_multianno.filter2.txt | cut -d: -f 2 | sed '1d' | awk -F ',' 'BEGIN{print "AF"}(\$1+\$2)!=0{print \$2/(\$1+\$2)}(\$1+\$2)==0{print "both is 0"}'| paste -d "\t"  ${libName}.anno.hg19_multianno.filter3.txt - > ${libName}.anno.hg19_multianno.final.txt
    perl ${params.filter_ch} ${libName}.anno.hg19_multianno.final.txt > ${libName}.anno.hg19_multianno.final.clinvar_hgmd.txt
    """
  }
}

process SV{
    tag "${libName}"
    cpus "${params.sv_cpu}"
    publishDir "${params.o}/04.SV1/${libName}/", mode: 'copy', overwrite: true
    input:
    set val(libName), file(bam), file(bai) from mapping_bam19
    output:
    set val(libName), file("*.sv.fusions.xls") into sv_result, sv_result1
    set val(libName), file("*.sv.fusions.pdf") into sv_pdf_result
    """
    ${params.svaba} run -t ${bam} -k ${Tbed500} -a ${libName} -G ${params.fa}
    grep -E "^#|ASDIS" ${libName}.svaba.sv.vcf > ${libName}.sv.tmp.vcf
    awk '\$1!~/#/' ${libName}.svaba.sv.vcf|awk '{\$2=\$2"\\t"\$2;OFS="\\t";print \$0}'|bedtools intersect -a - -b /fuer2/01.Pipeline/03.pipeline_test/whhe/01.ST_DNA/db/PON/fusion_gene.bed -wa |sort -u -V|awk '\$7!=0{print \$1"\\t"\$2"\\t"\$4"\\t"\$5"\\t"\$6"\\t"\$7"\\t"\$8"\\t"\$9"\\t"\$10"\\t"\$11}'|cat ${libName}.sv.tmp.vcf - > ${libName}.sv.tmp1.vcf
    perl ${params.sv2vcf} ${libName}.sv.tmp1.vcf|sort -u -V  > ${libName}.sv.tmp.avinput
    ${params.anno_var} -buildver hg19 ${libName}.sv.tmp.avinput ${params.anno_db}
    perl ${params.get_sv} ${libName}.sv.tmp.avinput.variant_function ${libName}.result > ${libName}.sv.xls
    num=`wc -l ${libName}.result|awk '{print \$1}'`
    if [[ \$num < "2" ]] ;then
       echo -e "#gene1\\tgene2\\tstrand1(gene/fusion)\\tstrand2(gene/fusion)\\tbreakpoint1\\tbreakpoint2\\tsite1\\tsite2\\ttype\\tdirection1\\tdirection2\\tsplit_reads1\\tsplit_reads2\\tdiscordant_mates\\tcoverage1\\tcoverage2\\tconfidence\\tclosest_genomic_breakpoint1\\tclosest_genomic_breakpoint2\\tfilters\\tfusion_transcript\\treading_frame\\tpeptide_sequence\\tread_identifiers" > ${libName}.sv.fusions.xls
       touch ${libName}.sv.fusions.pdf
       exit
    else
    awk -F '\\t' 'NR==FNR{a[\$5]=\$6;next}NR!=FNR&&FNR>1{print \$1"\\t"\$2"\\t"a[\$1]"/"a[\$1]"\\t"a[\$2]"/"a[\$2]"\\t"\$3"\\t"\$5"\\t"\$4"\\t"\$6"\\t-\\t"\$12"\\t"\$13"\\t"\$10"\\t0\\t0\\t"\$9"\\t"\$9"\\thigh\\t.\\t.\\t.\\t.\\tin-frame\\t"\$NF"\\t.\\t"\$1"\\t"\$2"\\tchr"\$3":"a[\$1]"\\tchr"\$5":"a[\$2]"\\t-\\t"\$10"\\t."}' /fuer2/Work/hewh/project/tumor_RNA/db/hg19/hg19_gene.info.txt ${libName}.result|sed '1i#gene1\\tgene2\\tstrand1(gene/fusion)\\tstrand2(gene/fusion)\\tbreakpoint1\\tbreakpoint2\\tsite1\\tsite2\\ttype\\tdirection1\\tdirection2\\tsplit_reads1\\tsplit_reads2\\tdiscordant_mates\\tcoverage1\\tcoverage2\\tconfidence\\tclosest_genomic_breakpoint1\\tclosest_genomic_breakpoint2\\tfiltersfusion_transcript\\treading_frame\\tpeptide_sequence\\tread_identifiers\\t#gene1\\tgene2\\tbreakpoint1:strand1(fusion)\\tbreakpoint2:strand2(fusion)\\ttype\\tTotal_reads\\tlocal_freq(%)' > ${libName}.fusion.tmp.txt
    sed -n '2,\$p' ${libName}.fusion.tmp.txt |cut -f5,6|sed 's/\\t/\\n/g'|awk -F ':' '{print \$1"\\t"\$2"\\t"\$2}'|bedtools intersect -a - -b /fuer2/02.Database/01.Database_project/hg19/hg19_cytoBand.txt -wao|cut -f7|awk 'BEGIN{print "Chromosomal band of breakpoint1\\tChromosomal band of breakpoint2"}{if(NR%2!=0)ORS="\\t";else ORS="\\n"}1' |awk '{split(\$1,a,".");split(\$2,b,".");print a[1]"\\t"b[1]}'|paste - ${libName}.fusion.tmp.txt |awk -F '\\t' '{split(\$7,a,":");split(\$8,b,":");\$7=\$7"\\t"a[1]\$1;\$8=\$8"\\t"b[1]\$2;for(i=3;i<=NF;i++){printf \$i"\\t"};print ""}'|sed 's/\\t\$//g' > ${libName}.fusion.tmp1.txt
    awk -F '\\t' 'NR==FNR{a[\$1"\\t"\$2]=1;next}NR!=FNR{if(FNR==1){print \$0"\\twhite_genelist"}else{if(\$1"\\t"\$2 in a){print \$0"\\t"1}else{print \$0"\\t"0}}}' /fuer2/01.Pipeline/03.pipeline_test/whhe/01.ST_DNA/db/PON/fusion_gene.new.txt ${libName}.fusion.tmp1.txt > ${libName}.fusion.tmp2.txt
    awk -F '\\t' 'NR>1{split(\$5,a,":");split(\$7,b,":");print a[1]"\\t"a[2]"\\t"a[2]"\\n"b[1]"\\t"b[2]"\\t"b[2]}' ${libName}.fusion.tmp2.txt|bedtools intersect -a - -b /fuer2/02.Database/01.Database_project/hg19/hg19.CDS_exon.bed -wao|awk '{if(\$5=="-1"){print \$1":"\$2"\\tout-frame"}else{print \$1":"\$2"\\tin-frame"}}'|awk -F '\\t' 'NR==FNR{a[\$1]=\$2;next}NR!=FNR{if(FNR==1){print \$0}else{\$24=a[\$7];OFS="\\t";print \$0}}' - ${libName}.fusion.tmp2.txt > ${libName}.sv.fusions.xls
    sh /fuer2/Work/hewh/project/tumor_RNA/plot_fusion/plot_fusion.sh ${libName}.sv.fusions.xls /fuer2/Work/hewh/project/tumor_RNA/plot_fusion/hg19.gtf
    mv fusions.pdf ${libName}.sv.fusions.pdf
    fi
    """
}

process fusion{
    tag "${libName}"
    cpus "${params.fusion_cpu}"
    publishDir "${params.o}/03.fusion1/${libName}/", mode: 'copy', overwrite: true
    input:
    set val(libName), file(bam) from mp_bam4
    output:
    set val(libName), file("*.fusions.xls") into fusion_out3, fusion_result
    set val(libName), file("*.fusions.pdf") into fusion_pdf_result
    """
    ${params.factera} -C -r 10 -o . ${libName}.markdup.sort.bam ${params.exons_bed} ${params.factera_fa} ${Tbed}
    num=`wc -l ${libName}.markdup.sort.factera.fusions.txt|awk '{print \$1}'`
    if [[ \$num < "1" ]] ;then
        echo -e "#gene1\\tgene2\\tstrand1(gene/fusion)\\tstrand2(gene/fusion)\\tbreakpoint1\\tbreakpoint2\\tsite1\\tsite2\\ttype\\tdirection1\\tdirection2\\tsplit_reads1\\tsplit_reads2\\tdiscordant_mates\\tcoverage1\\tcoverage2\\tconfidence\\tclosest_genomic_breakpoint1\\tclosest_genomic_breakpoint2\\tfilters\\tfusion_transcript\\treading_frame\\tpeptide_sequence\\tread_identifiers" > ${libName}.fusions.xls
        touch ${libName}.fusions.pdf
        exit
    else
    sed -n '2,\$p' ${libName}.markdup.sort.factera.fusions.txt|awk -F '\\t' 'NR==FNR{a[\$5]=\$6;next}NR!=FNR{if(a[\$2]~/+/){d1="downstream"}else{d1="upstream"};if(a[\$3]~/+/){d2="upstream"}else{d2="downstream"};dis=\$12*2-\$6-\$7;if(dis <1){dis=0};if(\$17<\$12*2){\$17=\$12*2};print \$2"\\t"\$3"\\t"a[\$2]"/"a[\$2]"\\t"a[\$3]"/"a[\$3]"\\t"\$4"\\t"\$5"\\tsplice-site\\tsplice-site\\t"\$1"\\t"d1"\\t"d2"\\t"\$6"\\t"\$7"\\t"dis"\\t"\$17"\\t"\$17"\\thigh\\t.\\t.\\t.\\t.\\tin-frame\\t"\$18"\\t."}' /fuer2/Work/hewh/project/tumor_RNA/db/hg19/hg19_gene.info.txt - |sed '1i#gene1\\tgene2\\tstrand1(gene/fusion)\\tstrand2(gene/fusion)\\tbreakpoint1\\tbreakpoint2\\tsite1\\tsite2\\ttype\\tdirection1\\tdirection2\\tsplit_reads1\\tsplit_reads2\\tdiscordant_mates\\tcoverage1\\tcoverage2\\tconfidence\\tclosest_genomic_breakpoint1\\tclosest_genomic_breakpoint2\\tfilters\\tfusion_transcript\\treading_frame\\tpeptide_sequence\\tread_identifiers' > ${libName}.fusion.plot.txt
    awk -F '\\t' 'NR==1{print \$0"\\t#gene1\\tgene2\\tbreakpoint1:strand1(fusion)\\tbreakpoint2:strand2(fusion)\\ttype\\tTotal_reads\\tlocal_freq(%)"}NR>1{split(\$3,a,"/");split(\$4,b,"/");total=\$12+\$13;gsub(/TRA/,"translocations",\$9);gsub(/INV/,"inversions",\$9);gsub(/DEL/,"deletions",\$9);OFS="\\t";print \$0"\\t"\$1"\\t"\$2"\\tchr"\$5":"a[1]"\\tchr"\$6":"b[1]"\\t"\$9"\\t"total"\\t."}' ${libName}.fusion.plot.txt >  ${libName}.fusion.tmp.txt
    sed -n '2,\$p' ${libName}.fusion.tmp.txt|cut -f5,6|sed 's/\\t/\\n/g'|awk -F ':' '{print \$1"\\t"\$2"\\t"\$2}'|bedtools intersect -a - -b /fuer2/02.Database/01.Database_project/hg19/hg19_cytoBand.txt -wao|cut -f7|awk 'BEGIN{print "Chromosomal band of breakpoint1\\tChromosomal band of breakpoint2"}{if(NR%2!=0)ORS="\\t";else ORS="\\n"}1' |awk '{split(\$1,a,".");split(\$2,b,".");print a[1]"\\t"b[1]}'|paste - ${libName}.fusion.tmp.txt|awk -F '\\t' '{split(\$7,a,":");split(\$8,b,":");\$7=\$7"\\t"a[1]\$1;\$8=\$8"\\t"b[1]\$2;for(i=3;i<=NF;i++){printf \$i"\\t"};print ""}'|sed 's/\\t\$//g' > ${libName}.fusion.tmp1.txt
   awk -F '\\t' 'NR==FNR{a[\$1"\\t"\$2]=1;next}NR!=FNR{if(FNR==1){print \$0"\\twhite_genelist"}else{if(\$1"\\t"\$2 in a){print \$0"\\t"1}else{print \$0"\\t"0}}}' ${params.Pon_fusion} ${libName}.fusion.tmp1.txt > ${libName}.fusion.tmp2.txt
    awk -F '\\t' 'NR>1{print \$5"\\n"\$7}' ${libName}.fusion.tmp2.txt|sed 's/:/\\t/g'|awk '{print \$1"\\t"\$2"\\t"\$2"\\t0\\t0"}' > ${libName}.avinput
    ${params.anno_var} -buildver hg19 ${libName}.avinput ${params.anno_db}
    awk -F '\\t' 'NR==FNR{a[\$3":"\$4]=\$1;next}NR!=FNR{if(FNR==1){print \$0}else{\$9=a[\$5];\$10=a[\$7];OFS="\\t";print \$0}}' ${libName}.avinput.variant_function ${libName}.fusion.tmp2.txt > ${libName}.fusion.tmp3.txt
    ${params.bedtools} intersect -a ${libName}.avinput -b /fuer2/02.Database/01.Database_project/hg19/hg19.CDS_exon.bed -wao | awk '{if(\$7=="-1"){print \$1":"\$2"\\tout-frame"}else{print \$1":"\$2"\\tin-frame"}}'|awk -F '\\t' 'NR==FNR{a[\$1]=\$2;next}NR!=FNR{if(FNR==1){print \$0}else{\$24=a[\$7];OFS="\\t";print \$0}}' - ${libName}.fusion.tmp3.txt > ${libName}.fusions.xls
    sh /fuer2/Work/hewh/project/tumor_RNA/plot_fusion/plot_fusion.sh ${libName}.fusion.plot.txt /fuer2/Work/hewh/project/tumor_RNA/plot_fusion/hg19.gtf
    mv fusions.pdf ${libName}.fusions.pdf
    fi
    """
  }

process CNV{
    tag "${libName}"
    cpus "${params.cnv_cpu}"
    input:
    set val(libName), file(bam), file(bai) from mapping_bam15
    output:
    set val(libName), file("*call.final.filter.cnv") into cnv_out2
    set val(libName), file("*.cnv_gene.info.xls") into cnv_result
    set val(libName), file("*.cnv_cover.pdf") into cnv_pdf_result
    """
    ${params.cnvkit} batch ${bam} -r ${params.cnn}  --drop-low-coverage --diagram --scatter --output-dir ./ -p ${params.cnv_cpu} 
    ${params.cnvkit} call ${libName}.markdup.sort.bintest.cns -o ${libName}.BH.call.cnv --purity 1
    ${params.cnvkit} call ${libName}.markdup.sort.call.cns -o ${libName}.call.cnv --purity 1
    awk -F '\\t' '{len=\$3-\$2+1;if(NR==1 || len>=1000 && \$6!=2 && \$8<=0.5){print \$0}}' ${libName}.call.cnv > ${libName}.call.filter.cnv
    awk -F '\\t' '{if(NR==1){print \$0"\\tcnv_num"}else{cn=sprintf("%.1f",2*2^\$5);print \$0"\\t"cn}}' ${libName}.call.filter.cnv > ${libName}.call.final.filter.cnv
    awk -F '\\t' 'NR==FNR{split(\$4,a,",");for(i in a){b[a[i]]=\$11}}NR!=FNR{if(\$1 in b){print \$1"\\t"b[\$1]}else{print \$1"\\t2.0"}}' ${libName}.call.final.filter.cnv /fuer2/01.Pipeline/03.pipeline_test/whhe/01.ST_DNA/db/PON/cnv_gene.txt|awk 'NR>1' |sed '1iGENE\tCN\tCNV type' > ${libName}.cnv_gene.info.xls
    sh ${params.gene_cnv} ${libName}.markdup.sort.cnr
    sed -n '2,\$p' ${libName}.markdup.sort.cnr|sort -V|awk -F '\\t' 'NR>1{\$6=2*2^\$6;if(\$6>6){\$6=6};gsub(/X/,"23",\$1);gsub(/Y/,"24",\$1);print \$4"\\t"\$1"\\t"\$2"\\t"\$6}'|sed '1iGENE\\tCHR\\tBP\\tCN' > cnv.plot.txt
    awk -F '\\t' 'NR==FNR{a[\$1]=\$2;next}NR!=FNR{if(\$1 in a){if(a[\$1] >=2.8){if(\$4>a[\$1]+1){\$4=a[\$1]+1};if(\$4<1.2){\$4=2}}if(a[\$1] <=1.2){if(\$4>2.8){\$4=2}if(\$4<0.1){\$4=0.1}}if(a[\$1]<2.8 && a[\$1]>1.2){if(\$4<1.2 || \$4>2.8){\$4=2}}}else{if(\$4>2.8 || \$4<1.2){\$4=2}};OFS="\\t";print \$0}' ${libName}.cnv_gene.info.xls cnv.plot.txt|sed 's/\\tBP\\t2/\\tBP\\tCN/' > cnv.plot.new.txt
    /fuer2/home/whhe/miniconda3/bin/Rscript ${params.cnv_plot} cnv.plot.new.txt ${params.Pon_gene} ${params.chr_center} ${params.chr_pos} ${libName}
    """
}

workflow.onComplete {
    if( workflow.success ) {
        File f = new File("./pipe.Done")
        f.write("Pipeline completed at: $workflow.complete"+"\n") 
        f.append("Execution status:      ${ workflow.success ? 'OK' : 'failed' }" + "\n")
        f.append("User:                  $workflow.userName" + "\n")
        f.append("Launch time:           ${workflow.start.format('yyyy-MMM-dd HH:mm:ss')}" + "\n")
        f.append("Ending time:           ${workflow.complete.format('yyyy-MMM-dd HH:mm:ss')}" + "\n")
        f.append("Duration:              $workflow.duration" + "\n")
        f.append("Total CPU-Hours:       ${workflow.stats.computeTimeFmt ?: '-'}" + "\n")
        f.append("Tasks stats:           Succeeded ${workflow.stats.succeedCountFmt}; Cached ${workflow.stats.cachedCountFmt}; Ignored ${workflow.stats.ignoredCountFmt}; Failed ${workflow.stats.failedCountFmt}"  + "\n")
    }else{
        File f = new File("./pipe.Failed")
        f.write("Pipeline failed at: $workflow.complete"+"\n")
        f.append("Execution status:   ${ workflow.success ? 'OK' : 'failed' }" + "\n")
        f.append("User:               $workflow.userName" + "\n")
        f.append("Launch time:        ${workflow.start.format('yyyy-MMM-dd HH:mm:ss')}" + "\n")
        f.append("Ending time:        ${workflow.complete.format('yyyy-MMM-dd HH:mm:ss')}" + "\n")
        f.append("Duration:           $workflow.duration" + "\n")
        f.append("Total CPU-Hours:    ${workflow.stats.computeTimeFmt ?: '-'}" + "\n")
        f.append("Tasks stats:        Succeeded ${workflow.stats.succeedCountFmt}; Cached ${workflow.stats.cachedCountFmt}; Ignored ${workflow.stats.ignoredCountFmt}; Failed ${workflow.stats.failedCountFmt}" + "\n")
        f.append("ERROR message:\n" + "  ${workflow.errorMessage}" + "\n")
        f.append("ERROR report:\n" + "  ${workflow.errorReport}" + "\n")
    }
}
