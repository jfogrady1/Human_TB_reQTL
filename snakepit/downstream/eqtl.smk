autosome = [str(i) for i in range(1, 23)]
#rule all:
     #input:
        #expand("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_{timepoint}.txt", timepoint = config["TIME"]),
        #expand("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_TPM_{timepoint}.txt", timepoint = config["TIME"]),
        #'/home/workspace/jogrady/heQTL/data/covariate/Genotypes.PCA_eigenvect.txt',
        #'/home/workspace/jogrady/heQTL/data/ref_genome/TSS_annotation.bed.gz',
        #expand('/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/TMM_INT_counts_{timepoint}.bed.gz', timepoint = config["TIME"]),
        #expand('/home/workspace/jogrady/heQTL/data/covariate/Scree_plot_PCA_{timepoint}.pdf', timepoint = config["TIME"]),
        #expand('/home/workspace/jogrady/heQTL/data/covariate/{timepoint}_eqtlcovs.txt', timepoint = config["TIME"]),
        #multiext("/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", ".bed", ".bim", ".fam"),          
        #expand("/home/workspace/jogrady/heQTL/results/eQTL/{timepoint}.{window}.cis_qtl.txt.gz", timepoint = config["TIME"], window = config["cis_window"]),
        #expand("/home/workspace/jogrady/heQTL/results/eQTL/{timepoint}.{window}.cis_qtl_fdr0.1.txt", window = config["cis_window"], timepoint = config["TIME"]),
        #expand('/home/workspace/jogrady/heQTL/results/eQTL/{timepoint}.{window}.cis_qtl_pairs.chr{chr}.parquet', timepoint = config["TIME"], window = config["cis_window"], chr = autosome),
        #expand('/home/workspace/jogrady/heQTL/results/eQTL/{timepoint}.{window}.cis_qtl_pairs.chr{chr}.txt.gz', timepoint = config["TIME"], window = config["cis_window"], chr = autosome),
        #expand("/home/workspace/jogrady/heQTL/results/eQTL/{timepoint}.50000.cis_qtl_pairs.ALL.txt", timepoint = config["TIME"]),
        #"/home/workspace/jogrady/heQTL/results/eQTL/numbers_stats_cis.txt"

rule genotype_PCA_covariate:
     input:
        vcf = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf.gz',
        script = '/home/workspace/jogrady/heQTL/scripts/Genotype_PCA_covariate.R'

     output:
        eigenvec = '/home/workspace/jogrady/heQTL/data/covariate/Genotypes.PCA_eigenvect.txt',
        eigenval = '/home/workspace/jogrady/heQTL/data/covariate/Genotypes.PCA_eigenval.txt',
        unzipped = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf'
     shell:
        '''
         bcftools view {input.vcf} -Ov -o {output.unzipped}
         Rscript {input.script} {output.unzipped}  /home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Genotypes.ccm.gds {output.eigenvec} {output.eigenval}
        '''
rule TPM:
     input: 
        counts = '/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_{timepoint}.txt', 
        annotation = '/home/workspace/jogrady/heQTL/data/ref_genome/gencode.v43.annotation.gtf',
        script = '/home/workspace/jogrady/heQTL/scripts/TPM_normalisation.R'
     output: 
        tpm = '/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_TPM_{timepoint}.txt'

     shell:
        '''
         Rscript {input.script} {input.counts} {input.annotation} {output.tpm}         
        '''

rule gtf2bed:
    input: 
        annotation = '/home/workspace/jogrady/heQTL/data/ref_genome/gencode.v43.annotation.gtf'
    output: 
        bed = '/home/workspace/jogrady/heQTL/data/ref_genome/TSS_annotation.bed.gz'

    shell:
       '''
       cat {input.annotation} | awk 'OFS="\t" {{if ($3=="gene") {{if ($7 == "+") {{print $1,$4-1,$4,$10}} else {{print $1,$5-1,$5,$10}}}}}}' | tr -d '";' | bgzip -c > {output.bed}
       '''

rule Filtering_TMM_Norm:
     input:
         counts = '/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_{timepoint}.txt',
         tpm = '/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_TPM_{timepoint}.txt',
         script = '/home/workspace/jogrady/heQTL/scripts/TMM_normalisation.R',
         tss = '/home/workspace/jogrady/heQTL/data/ref_genome/TSS_annotation.bed.gz',
         geno = '/home/workspace/jogrady/heQTL/data/covariate/Genotypes.PCA_eigenvect.txt',
         cov = '/home/workspace/jogrady/heQTL/data/covariate/{timepoint}_metadata.txt'
     output:
         bed = '/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/TMM_INT_counts_{timepoint}.bed.gz',
         scree = '/home/workspace/jogrady/heQTL/data/covariate/Scree_plot_PCA_{timepoint}.pdf',
         covariate = '/home/workspace/jogrady/heQTL/data/covariate/{timepoint}_eqtlcovs.txt',
     resources:
        mem_mb = 6000,
        threads = 10
     shell:
       '''
       Rscript {input.script} {input.counts} {input.tpm} {input.tss} /home/workspace/jogrady/heQTL/work/RNA_seq/quantification/TMM_INT_counts_{wildcards.timepoint}.bed {output.scree} {input.cov} {input.geno} {output.covariate}
       bgzip /home/workspace/jogrady/heQTL/work/RNA_seq/quantification/TMM_INT_counts_{wildcards.timepoint}.bed
       tabix -p bed /home/workspace/jogrady/heQTL/work/RNA_seq/quantification/TMM_INT_counts_{wildcards.timepoint}.bed.gz
       '''
rule prepareplink:
    input:
        vcf = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf'
    output:
        plink = multiext("/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", ".bed", ".bim", ".fam")
    
    params:
       prefix = "/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/eQTL_genotypes_R0.6"
    shell:
        '''
        plink --vcf {input.vcf} --output-chr chrM --double-id --keep-allele-order --autosome --make-bed --out {params.prefix}
        '''

rule qtl_mapping_permute:
    input:
        expression_bed='/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/TMM_INT_counts_{timepoint}.bed.gz', # Phenotypes
        covariates_file='/home/workspace/jogrady/heQTL/data/covariate/{timepoint}_eqtlcovs.txt', # Covariates
        #pyscript = '/home/workspace/jogrady/heQTL/scripts/parquet2txt.py'
    output:
        outputdir="/home/workspace/jogrady/heQTL/results/eQTL/{timepoint}.{window}.cis_qtl.txt.gz"
    resources:
        mem_mb = 6000,
        threads = 40
    params:
       plink_prefix_path="/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", # Genotypes
       prefix = "/home/workspace/jogrady/heQTL/results/eQTL/{timepoint}.{window}",
       cis_window_size = "{window}"
    shell:
        '''
         python3 -m tensorqtl {params.plink_prefix_path} {input.expression_bed} {params.prefix} \
            --mode cis \
            --window {params.cis_window_size} \
            --seed 1856 \
            --permutations 10000 \
            --covariates {input.covariates_file} \
            --fdr 0.1
        '''
          
rule eGenedetection:
    input:
        eqtl_permute='/home/workspace/jogrady/heQTL/results/eQTL/{timepoint}.{window}.cis_qtl.txt.gz', # Phenotypes
        script='/home/workspace/jogrady/heQTL/scripts/eGene_detection.R'
    output:
        eGenes= '/home/workspace/jogrady/heQTL/results/eQTL/{timepoint}.{window}.cis_qtl_fdr0.1.txt'
    shell:
        '''
        Rscript {input.script} {input.eqtl_permute} {output.eGenes} 0.1
        '''

rule nominal_mappinngT0:
    input:
        expression_bed='/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/TMM_INT_counts_T0.bed.gz', # Phenotypes
        covariates_file='/home/workspace/jogrady/heQTL/data/covariate/T0_eqtlcovs.txt', # Covariates
    output:
        parquet_files= '/home/workspace/jogrady/heQTL/results/eQTL/T0.{window}.cis_qtl_pairs.chr{chr}.parquet'
    resources:
        mem_mb = 6000,
        threads = 40
    params:
       plink_prefix_path="/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", # Genotypes
       prefix = "/home/workspace/jogrady/heQTL/results/eQTL/T0.{window}",
       cis_window_size = "{window}"
    shell:
        '''
        python3 -m tensorqtl {params.plink_prefix_path} {input.expression_bed} {params.prefix} \
            --covariates {input.covariates_file} \
            --window {params.cis_window_size} \
            --mode cis_nominal
        '''
rule nominal_mappinngT1:
    input:
        expression_bed='/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/TMM_INT_counts_T1.bed.gz', # Phenotypes
        covariates_file='/home/workspace/jogrady/heQTL/data/covariate/T1_eqtlcovs.txt', # Covariates
    output:
        parquet_files= '/home/workspace/jogrady/heQTL/results/eQTL/T1.{window}.cis_qtl_pairs.chr{chr}.parquet'
    resources:
        mem_mb = 6000,
        threads = 40
    params:
       plink_prefix_path="/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", # Genotypes
       prefix = "/home/workspace/jogrady/heQTL/results/eQTL/T1.{window}",
       cis_window_size = "{window}"
    shell:
        '''
        python3 -m tensorqtl {params.plink_prefix_path} {input.expression_bed} {params.prefix} \
            --covariates {input.covariates_file} \
            --window {params.cis_window_size} \
            --mode cis_nominal
        '''
rule nominal_mappinngT2:
    input:
        expression_bed='/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/TMM_INT_counts_T2.bed.gz', # Phenotypes
        covariates_file='/home/workspace/jogrady/heQTL/data/covariate/T2_eqtlcovs.txt', # Covariates
    output:
        parquet_files= '/home/workspace/jogrady/heQTL/results/eQTL/T2.{window}.cis_qtl_pairs.chr{chr}.parquet'
    resources:
        mem_mb = 6000,
        threads = 40
    params:
       plink_prefix_path="/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", # Genotypes
       prefix = "/home/workspace/jogrady/heQTL/results/eQTL/T2.{window}",
       cis_window_size = "{window}"
    shell:
        '''
        python3 -m tensorqtl {params.plink_prefix_path} {input.expression_bed} {params.prefix} \
            --covariates {input.covariates_file} \
            --window {params.cis_window_size} \
            --mode cis_nominal
        '''
rule nominal_mappinngT3:
    input:
        expression_bed='/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/TMM_INT_counts_T3.bed.gz', # Phenotypes
        covariates_file='/home/workspace/jogrady/heQTL/data/covariate/T3_eqtlcovs.txt', # Covariates,
    output:
        parquet_files= '/home/workspace/jogrady/heQTL/results/eQTL/T3.{window}.cis_qtl_pairs.chr{chr}.parquet'
    resources:
        mem_mb = 6000,
        threads = 40
    params:
       plink_prefix_path="/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", # Genotypes
       prefix = "/home/workspace/jogrady/heQTL/results/eQTL/T3.{window}",
       cis_window_size = "{window}"
    shell:
        '''
        python3 -m tensorqtl {params.plink_prefix_path} {input.expression_bed} {params.prefix} \
            --covariates {input.covariates_file} \
            --window {params.cis_window_size} \
            --mode cis_nominal
        '''

rule nominal_mappinngT4:
    input:
        expression_bed='/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/TMM_INT_counts_T4.bed.gz', # Phenotypes
        covariates_file='/home/workspace/jogrady/heQTL/data/covariate/T4_eqtlcovs.txt', # Covariates
    output:
        parquet_files= '/home/workspace/jogrady/heQTL/results/eQTL/T4.{window}.cis_qtl_pairs.chr{chr}.parquet'
    resources:
        mem_mb = 6000,
        threads = 40
    params:
       plink_prefix_path="/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", # Genotypes
       prefix = "/home/workspace/jogrady/heQTL/results/eQTL/T4.{window}",
       cis_window_size = "{window}"
    shell:
        '''
        python3 -m tensorqtl {params.plink_prefix_path} {input.expression_bed} {params.prefix} \
            --covariates {input.covariates_file} \
            --window {params.cis_window_size} \
            --mode cis_nominal
        '''



rule parquet2py:
    input:
        pyscript = '/home/workspace/jogrady/heQTL/scripts/parquet2txt.py',
        parquet = '/home/workspace/jogrady/heQTL/results/eQTL/{timepoint}.50000.cis_qtl_pairs.chr{chr}.parquet'

    output:
        txt = '/home/workspace/jogrady/heQTL/results/eQTL/{timepoint}.50000.cis_qtl_pairs.chr{chr}.txt.gz'

    shell:
        '''     
        python3 {input.pyscript} {input.parquet} {output.txt}
        '''

rule merged_nominal:
    input:
        T0_nominal = expand("/home/workspace/jogrady/heQTL/results/eQTL/T0.50000.cis_qtl_pairs.chr{autosome}.txt.gz", autosome = autosome),
        T1_nominal = expand("/home/workspace/jogrady/heQTL/results/eQTL/T1.50000.cis_qtl_pairs.chr{autosome}.txt.gz", autosome = autosome),
        T2_nominal = expand("/home/workspace/jogrady/heQTL/results/eQTL/T2.50000.cis_qtl_pairs.chr{autosome}.txt.gz", autosome = autosome),
        T3_nominal = expand("/home/workspace/jogrady/heQTL/results/eQTL/T3.50000.cis_qtl_pairs.chr{autosome}.txt.gz", autosome = autosome),
        T4_nominal = expand("/home/workspace/jogrady/heQTL/results/eQTL/T4.50000.cis_qtl_pairs.chr{autosome}.txt.gz", autosome = autosome)

    output:
        T0_nominal_merged = "/home/workspace/jogrady/heQTL/results/eQTL/T0.50000.cis_qtl_pairs.ALL.txt.gz",
        T1_nominal_merged = "/home/workspace/jogrady/heQTL/results/eQTL/T1.50000.cis_qtl_pairs.ALL.txt.gz",
        T2_nominal_merged = "/home/workspace/jogrady/heQTL/results/eQTL/T2.50000.cis_qtl_pairs.ALL.txt.gz",
        T3_nominal_merged = "/home/workspace/jogrady/heQTL/results/eQTL/T3.50000.cis_qtl_pairs.ALL.txt.gz",
        T4_nominal_merged = "/home/workspace/jogrady/heQTL/results/eQTL/T4.50000.cis_qtl_pairs.ALL.txt.gz"

    shell:
        '''
        cat {input.T0_nominal} > {output.T0_nominal_merged}
        cat {input.T1_nominal} > {output.T1_nominal_merged}
        cat {input.T2_nominal} > {output.T2_nominal_merged}
        cat {input.T3_nominal} > {output.T3_nominal_merged}
        cat {input.T4_nominal} > {output.T4_nominal_merged}
        '''


rule merged_unzip:
    input:
        T0_nominal_merged = "/home/workspace/jogrady/heQTL/results/eQTL/T0.50000.cis_qtl_pairs.ALL.txt.gz",
        T1_nominal_merged = "/home/workspace/jogrady/heQTL/results/eQTL/T1.50000.cis_qtl_pairs.ALL.txt.gz",
        T2_nominal_merged = "/home/workspace/jogrady/heQTL/results/eQTL/T2.50000.cis_qtl_pairs.ALL.txt.gz",
        T3_nominal_merged = "/home/workspace/jogrady/heQTL/results/eQTL/T3.50000.cis_qtl_pairs.ALL.txt.gz",
        T4_nominal_merged = "/home/workspace/jogrady/heQTL/results/eQTL/T4.50000.cis_qtl_pairs.ALL.txt.gz"

    output:
        T0_unzipped_merged = "/home/workspace/jogrady/heQTL/results/eQTL/T0.50000.cis_qtl_pairs.ALL.txt",
        T1_unzipped_merged = "/home/workspace/jogrady/heQTL/results/eQTL/T1.50000.cis_qtl_pairs.ALL.txt",
        T2_unzipped_merged = "/home/workspace/jogrady/heQTL/results/eQTL/T2.50000.cis_qtl_pairs.ALL.txt",
        T3_unzipped_merged = "/home/workspace/jogrady/heQTL/results/eQTL/T3.50000.cis_qtl_pairs.ALL.txt",
        T4_unzipped_merged = "/home/workspace/jogrady/heQTL/results/eQTL/T4.50000.cis_qtl_pairs.ALL.txt"

    shell:
        '''
        gunzip -c {input.T0_nominal_merged} > {output.T0_unzipped_merged}
        gunzip -c {input.T1_nominal_merged} > {output.T1_unzipped_merged}
        gunzip -c {input.T2_nominal_merged} > {output.T2_unzipped_merged}
        gunzip -c {input.T3_nominal_merged} > {output.T3_unzipped_merged}
        gunzip -c {input.T4_nominal_merged} > {output.T4_unzipped_merged}
        '''

rule number_eQTLs:
    input:
        permute = expand("/home/workspace/jogrady/heQTL/results/eQTL/{timepoint}.50000.cis_qtl_fdr0.1.txt", timepoint = config["TIME"]),
        nominal = expand("/home/workspace/jogrady/heQTL/results/eQTL/{timepoint}.50000.cis_qtl_pairs.ALL.txt", timepoint = config["TIME"]),
        script = "/home/workspace/jogrady/heQTL/scripts/cis-eQTL_numbers.R"
    output:
        threshold_violin = "/home/workspace/jogrady/heQTL/results/eQTL/results/violin_plot_pval_thresholds.pdf",
        eGene_tested = "/home/workspace/jogrady/heQTL/results/eQTL/results/upset_cis_egenes_tested.pdf",
        upset_plot = "/home/workspace/jogrady/heQTL/results/eQTL/results/upset_cis_egenes_fdr_0.1.pdf",
        stats = "/home/workspace/jogrady/heQTL/results/eQTL/numbers_stats_cis.txt",
        TWAS_qtls = expand("/home/workspace/jogrady/heQTL/results/eQTL/{timepoint}_cis_qtl_pairs_TWAS_input.txt", timepoint = config["TIME"])
    shell:
        '''
        Rscript {input.script} {input.permute} {input.nominal} {output.threshold_violin} {output.eGene_tested} {output.stats} {output.TWAS_qtls} {output.upset_plot}
        '''

