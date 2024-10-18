# Main Snakefile
configfile: "/home/workspace/jogrady/heQTL/config/deepvariant.yaml"  # Primary config file (used for T0 - T4)

import yaml
autosomes = [str(i) for i in range(1, 23)]
ids = [1,5,6,7]
# Load the second config file for the merged BAMs
with open("/home/workspace/jogrady/heQTL/config/Alignment.yaml") as f:
    config2 = yaml.safe_load(f)  # This contains keys for merged BAMs



include: 'snakepit/variant_call/T0_deepvariant.smk'
include: 'snakepit/variant_call/T1_deepvariant.smk'
include: 'snakepit/variant_call/T2_deepvariant.smk'
include: 'snakepit/variant_call/T3_deepvariant.smk'
include: 'snakepit/variant_call/T4_deepvariant.smk'
include: 'snakepit/variant_call/Concordance.smk'
include: 'snakepit/variant_call/mergebam.smk'
include: 'snakepit/variant_call/mosdepth_coverage.smk'
include: 'snakepit/variant_call/Merged_deepvariant.smk'
include: 'snakepit/variant_call/giab_comparison.smk'
include: 'snakepit/variant_call/imputation.smk'
include: 'snakepit/variant_call/imputation_QC.smk'
include: 'snakepit/variant_call/sample_mismatch_eQTL.smk'
include: 'snakepit/variant_call/WASP_mapping_filtering.smk'



rule all:
    input:
        #T0 - T4 - note these are all just single files
        expand('/home/workspace/jogrady/heQTL/work/RNA_seq/variant_call/T0/{region}.{preset}.vcf.gz',region=config['regions'],preset=config['GL_config']),
        expand('/home/workspace/jogrady/heQTL/work/RNA_seq/variant_call/T1/{region}.{preset}.vcf.gz',region=config['regions'],preset=config['GL_config']),
        expand('/home/workspace/jogrady/heQTL/work/RNA_seq/variant_call/T2/{region}.{preset}.vcf.gz',region=config['regions'],preset=config['GL_config']),
        expand('/home/workspace/jogrady/heQTL/work/RNA_seq/variant_call/T3/{region}.{preset}.vcf.gz',region=config['regions'],preset=config['GL_config']),
        expand('/home/workspace/jogrady/heQTL/work/RNA_seq/variant_call/T4/{region}.{preset}.vcf.gz',region=config['regions'],preset=config['GL_config']),
        # Concordance
        expand('/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/{timepoint}_Genotypes.GT.FORMAT', timepoint=config2["TIME"]),
        '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T0_V_T1.pdf',
        '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T1_V_T0.pdf',
        '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T2_V_T0.pdf',
        '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T3_V_T0.pdf',
        '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T4_V_T0.pdf',
        # merged bams
        expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/{id}.bam.bai', id = config2["merged_patient_IDS"]),
        expand('/home/workspace/jogrady/heQTL/work/RNA_seq/variant_call/mergedcall/deepvariant/{sample}.{region}.g.vcf.gz',sample=config["samples_merged"], region = config["regions"], allow_missing=True),
        # coverage check (on merged files)
        expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/{id}_coverage.mosdepth.region.dist.txt', id=config2["merged_patient_IDS"]),
        expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/coding_regions.html'),
        expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/{id}_coverage_WG.mosdepth.global.dist.txt', id=config2["merged_patient_IDS"]),
        expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/WG_regions.html'),
        # Coverage
        expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/{id}_coverage.mosdepth.region.dist.txt', id=config2["merged_patient_IDS"]),
        expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/coding_regions.html'),
        expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/{id}_coverage_WG.mosdepth.global.dist.txt', id=config2["merged_patient_IDS"]),
        expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/WG_regions.html'),
        # Final variant call (on merged files)
        expand('/home/workspace/jogrady/heQTL/work/RNA_seq/variant_call/mergedcall/{region}.{preset}.vcf.gz',region=config['regions'],preset=config['GL_config']),
        # GIAB
        expand('/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/variantcomparison/HG00{id}_calls.txt', id = ids),
        '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/variantcomparison/DV_nerged_raw_calls.txt',
        '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/variantcomparison/DV_raw_GIAB_concordance.pdf',
        # Imputation
        '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/all.filtered.vcf.gz',
        expand('/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/DV_chr{chromosome}_IMPUTED.vcf.gz.tbi', chromosome=autosomes),
        expand('/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/all.filtered.chr{chromosomes}.bcf', chromosomes=autosomes),
        expand('/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/chr_{chromosome}_phased.bcf', chromosome=autosomes),
        expand('/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/DV_chr{chromosome}_IMPUTED.vcf.gz', chromosome=autosomes),
        '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_raw_IMPUTED.vcf.gz',
        '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf.gz',
        '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_IMPUTED_R0.3_filtered.vcf.gz',
        # Imputation_QC
        '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/chr14_karyogram_raw.pdf',
        '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/PCA_pruned_imputed_variants.pdf',
        '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Chr14_CDS.bed',
        '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Chr14.vcf',
        '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Typed_variants_chr14.bed',
        '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Imputed_variants_chr14.bed',
        '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Imputed_variants_R0.3_chr14.bed',
        expand('/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Overlap_Chr14_CDS_R0.{R}.txt', R = [6,3]),
        '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/SNP_density_all_chr.pdf',
        '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/chr14_karyogram_raw.pdf',
        # MBV
        expand("/home/workspace/jogrady/heQTL/results/eQTL/sample_check/{sample}.mbv_output.txt", sample = config2["IDS_0"]),
        expand("/home/workspace/jogrady/heQTL/results/eQTL/sample_check/{sample}.mbv_output.txt", sample = config2["IDS_1"]),
        expand("/home/workspace/jogrady/heQTL/results/eQTL/sample_check/{sample}.mbv_output.txt", sample = config2["IDS_2"]),
        expand("/home/workspace/jogrady/heQTL/results/eQTL/sample_check/{sample}.mbv_output.txt", sample = config2["IDS_3"]),
        expand("/home/workspace/jogrady/heQTL/results/eQTL/sample_check/{sample}.mbv_output.txt", sample = config2["IDS_4"]),
        #WASP (one sample only)
        "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/DV_snp_tab.h5",
        "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/intersecting_snps/SRR12609781.keep.bam",
        '/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/map2/SRR12609781.sort.bam',
        "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/filter_remapped/SRR12609781.keep.bam",
        "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/merge/SRR12609781.keep.merge.sort.bam",
        "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/dup_filtered/sorted/SRR12609781.keep.merge.sort.no_dups.bam",
        "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/dup_filtered/sorted/SRR12609781.mbv_output.txt"
