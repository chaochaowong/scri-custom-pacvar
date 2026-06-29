include { VEP_EXTRACT_CLINVAR_COSMIC } from '../../../modules/local/vep/extract_clinvar_cosmic'

workflow VEP_POST_ANNOTATION {
    take:
    ch_vcf_tbi

    main:
    VEP_EXTRACT_CLINVAR_COSMIC(ch_vcf_tbi)

    emit:
    clinvar = VEP_EXTRACT_CLINVAR_COSMIC.out.clinvar
    cosmic  = VEP_EXTRACT_CLINVAR_COSMIC.out.cosmic
}
