process VEP_EXTRACT_CLINVAR_COSMIC {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'python:3.12-slim'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${prefix}.clinvar_pathogenic.csv"), emit: clinvar
    tuple val(meta), path("${prefix}.cosmic.csv")            , emit: cosmic
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.file_name ?: meta.id}.vep"
    def args = task.ext.args ?: ''
    """
    extract_vep_clinvar_cosmic_csv.py \\
        ${vcf} \\
        --clinvar-output ${prefix}.clinvar_pathogenic.csv \\
        --cosmic-output ${prefix}.cosmic.csv \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.file_name ?: meta.id}.vep"
    """
    printf 'CHROM,POS,ID,REF,ALT,ALT_INDEX,QUAL,FILTER,HGVSg,Allele,Consequence,IMPACT,SYMBOL,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,HGVSc,HGVSp,Existing_variation,CLIN_SIG,SOMATIC,PHENO,PUBMED,CANONICAL,MANE_SELECT,AF,gnomADe_AF,gnomADg_AF,MAX_AF,MAX_AF_POPS,SAMPLE,GT,GQ,DP,AD,AD_REF,AD_ALT,REF_ALLELE_FRACTION,ALT_ALLELE_FRACTION,VAF,PS\\n' > ${prefix}.clinvar_pathogenic.csv
    printf 'CHROM,POS,ID,REF,ALT,ALT_INDEX,QUAL,FILTER,HGVSg,Allele,Consequence,IMPACT,SYMBOL,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,HGVSc,HGVSp,Existing_variation,CLIN_SIG,SOMATIC,PHENO,PUBMED,CANONICAL,MANE_SELECT,AF,gnomADe_AF,gnomADg_AF,MAX_AF,MAX_AF_POPS,SAMPLE,GT,GQ,DP,AD,AD_REF,AD_ALT,REF_ALLELE_FRACTION,ALT_ALLELE_FRACTION,VAF,PS\\n' > ${prefix}.cosmic.csv
    """
}
