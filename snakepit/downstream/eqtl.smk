autosome = [str(i) for i in range(1, 23)]
#rule all:
     #input:
        #expand("work/RNA_seq/quantification/count_matrix_ordered_{timepoint}.txt", timepoint = config["TIME"]),
        #expand("work/RNA_seq/quantification/count_matrix_TPM_{timepoint}.txt", timepoint = config["TIME"]),
        #'data/covariate/Genotypes.PCA_eigenvect.txt',
        #'data/ref_genome/TSS_annotation.bed.gz',
        #expand('work/RNA_seq/quantification/TMM_INT_counts_{timepoint}.bed.gz', timepoint = config["TIME"]),
        #expand('data/covariate/Scree_plot_PCA_{timepoint}.pdf', timepoint = config["TIME"]),
        #expand('data/covariate/{timepoint}_eqtlcovs.txt', timepoint = config["TIME"]),
        #multiext("work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", ".bed", ".bim", ".fam"),          
        #expand("results/eQTL/{timepoint}.{window}.cis_qtl.txt.gz", timepoint = config["TIME"], window = config["cis_window"]),
        #expand("results/eQTL/{timepoint}.{window}.cis_qtl_fdr0.1.txt", window = config["cis_window"], timepoint = config["TIME"]),
        #expand('results/eQTL/{timepoint}.{window}.cis_qtl_pairs.chr{chr}.parquet', timepoint = config["TIME"], window = config["cis_window"], chr = autosome),
        #expand('results/eQTL/{timepoint}.{window}.cis_qtl_pairs.chr{chr}.txt.gz', timepoint = config["TIME"], window = config["cis_window"], chr = autosome),
        #expand("results/eQTL/{timepoint}.50000.cis_qtl_pairs.ALL.txt", timepoint = config["TIME"]),
        #"results/eQTL/numbers_stats_cis.txt"

rule genotype_PCA_covariate:
     input:
        vcf = 'work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf.gz',
        script = 'scripts/Genotype_PCA_covariate.R'

     output:
        eigenvec = 'data/covariate/Genotypes.PCA_eigenvect.txt',
        eigenval = 'data/covariate/Genotypes.PCA_eigenval.txt',
        unzipped = 'work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf'
     shell:
        '''
         bcftools view {input.vcf} -Ov -o {output.unzipped}
         Rscript {input.script} {output.unzipped}  work/DNA_seq/imputation/final/Genotypes.ccm.gds {output.eigenvec} {output.eigenval}
        '''
rule TPM:
     input: 
        counts = 'work/RNA_seq/quantification/count_matrix_ordered_{timepoint}.txt', 
        annotation = 'data/ref_genome/gencode.v43.annotation.gtf',
        script = 'scripts/TPM_normalisation.R'
     output: 
        tpm = 'work/RNA_seq/quantification/count_matrix_TPM_{timepoint}.txt'

     shell:
        '''
         Rscript {input.script} {input.counts} {input.annotation} {output.tpm}         
        '''

rule gtf2bed:
    input: 
        annotation = 'data/ref_genome/gencode.v43.annotation.gtf'
    output: 
        bed = 'data/ref_genome/TSS_annotation.bed.gz'

    shell:
       '''
       cat {input.annotation} | awk 'OFS="\t" {{if ($3=="gene") {{if ($7 == "+") {{print $1,$4-1,$4,$10}} else {{print $1,$5-1,$5,$10}}}}}}' | tr -d '";' | bgzip -c > {output.bed}
       '''

rule Filtering_TMM_Norm:
     input:
         counts = 'work/RNA_seq/quantification/count_matrix_ordered_{timepoint}.txt',
         tpm = 'work/RNA_seq/quantification/count_matrix_TPM_{timepoint}.txt',
         script = 'scripts/TMM_normalisation.R',
         tss = 'data/ref_genome/TSS_annotation.bed.gz',
         geno = 'data/covariate/Genotypes.PCA_eigenvect.txt',
         cov = 'data/covariate/{timepoint}_metadata.txt'
     output:
         bed = 'work/RNA_seq/quantification/TMM_INT_counts_{timepoint}.bed.gz',
         scree = 'data/covariate/Scree_plot_PCA_{timepoint}.pdf',
         covariate = 'data/covariate/{timepoint}_eqtlcovs.txt',
     resources:
        mem_mb = 6000,
        threads = 10
     shell:
       '''
       Rscript {input.script} {input.counts} {input.tpm} {input.tss} work/RNA_seq/quantification/TMM_INT_counts_{wildcards.timepoint}.bed {output.scree} {input.cov} {input.geno} {output.covariate}
       bgzip work/RNA_seq/quantification/TMM_INT_counts_{wildcards.timepoint}.bed
       tabix -p bed work/RNA_seq/quantification/TMM_INT_counts_{wildcards.timepoint}.bed.gz
       '''
rule prepareplink:
    input:
        vcf = 'work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf'
    output:
        plink = multiext("work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", ".bed", ".bim", ".fam")
    
    params:
       prefix = "work/DNA_seq/imputation/final/eQTL_genotypes_R0.6"
    shell:
        '''
        plink --vcf {input.vcf} --output-chr chrM --double-id --keep-allele-order --autosome --make-bed --out {params.prefix}
        '''

rule qtl_mapping_permute:
    input:
        expression_bed='work/RNA_seq/quantification/TMM_INT_counts_{timepoint}.bed.gz', # Phenotypes
        covariates_file='data/covariate/{timepoint}_eqtlcovs.txt', # Covariates
        #pyscript = 'scripts/parquet2txt.py'
    output:
        outputdir="results/eQTL/{timepoint}.{window}.cis_qtl.txt.gz"
    resources:
        mem_mb = 6000,
        threads = 40
    params:
       plink_prefix_path="work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", # Genotypes
       prefix = "results/eQTL/{timepoint}.{window}",
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
        eqtl_permute='results/eQTL/{timepoint}.{window}.cis_qtl.txt.gz', # Phenotypes
        script='scripts/eGene_detection.R'
    output:
        eGenes= 'results/eQTL/{timepoint}.{window}.cis_qtl_fdr0.1.txt'
    shell:
        '''
        Rscript {input.script} {input.eqtl_permute} {output.eGenes} 0.1
        '''

rule nominal_mappinngT0:
    input:
        expression_bed='work/RNA_seq/quantification/TMM_INT_counts_T0.bed.gz', # Phenotypes
        covariates_file='data/covariate/T0_eqtlcovs.txt', # Covariates
    output:
        parquet_files= 'results/eQTL/T0.{window}.cis_qtl_pairs.chr{chr}.parquet'
    resources:
        mem_mb = 6000,
        threads = 40
    params:
       plink_prefix_path="work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", # Genotypes
       prefix = "results/eQTL/T0.{window}",
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
        expression_bed='work/RNA_seq/quantification/TMM_INT_counts_T1.bed.gz', # Phenotypes
        covariates_file='data/covariate/T1_eqtlcovs.txt', # Covariates
    output:
        parquet_files= 'results/eQTL/T1.{window}.cis_qtl_pairs.chr{chr}.parquet'
    resources:
        mem_mb = 6000,
        threads = 40
    params:
       plink_prefix_path="work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", # Genotypes
       prefix = "results/eQTL/T1.{window}",
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
        expression_bed='work/RNA_seq/quantification/TMM_INT_counts_T2.bed.gz', # Phenotypes
        covariates_file='data/covariate/T2_eqtlcovs.txt', # Covariates
    output:
        parquet_files= 'results/eQTL/T2.{window}.cis_qtl_pairs.chr{chr}.parquet'
    resources:
        mem_mb = 6000,
        threads = 40
    params:
       plink_prefix_path="work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", # Genotypes
       prefix = "results/eQTL/T2.{window}",
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
        expression_bed='work/RNA_seq/quantification/TMM_INT_counts_T3.bed.gz', # Phenotypes
        covariates_file='data/covariate/T3_eqtlcovs.txt', # Covariates,
    output:
        parquet_files= 'results/eQTL/T3.{window}.cis_qtl_pairs.chr{chr}.parquet'
    resources:
        mem_mb = 6000,
        threads = 40
    params:
       plink_prefix_path="work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", # Genotypes
       prefix = "results/eQTL/T3.{window}",
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
        expression_bed='work/RNA_seq/quantification/TMM_INT_counts_T4.bed.gz', # Phenotypes
        covariates_file='data/covariate/T4_eqtlcovs.txt', # Covariates
    output:
        parquet_files= 'results/eQTL/T4.{window}.cis_qtl_pairs.chr{chr}.parquet'
    resources:
        mem_mb = 6000,
        threads = 40
    params:
       plink_prefix_path="work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", # Genotypes
       prefix = "results/eQTL/T4.{window}",
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
        pyscript = 'scripts/parquet2txt.py',
        parquet = 'results/eQTL/{timepoint}.50000.cis_qtl_pairs.chr{chr}.parquet'

    output:
        txt = 'results/eQTL/{timepoint}.50000.cis_qtl_pairs.chr{chr}.txt.gz'

    shell:
        '''     
        python3 {input.pyscript} {input.parquet} {output.txt}
        '''

rule merged_nominal:
    input:
        T0_nominal = expand("results/eQTL/T0.50000.cis_qtl_pairs.chr{autosome}.txt.gz", autosome = autosome),
        T1_nominal = expand("results/eQTL/T1.50000.cis_qtl_pairs.chr{autosome}.txt.gz", autosome = autosome),
        T2_nominal = expand("results/eQTL/T2.50000.cis_qtl_pairs.chr{autosome}.txt.gz", autosome = autosome),
        T3_nominal = expand("results/eQTL/T3.50000.cis_qtl_pairs.chr{autosome}.txt.gz", autosome = autosome),
        T4_nominal = expand("results/eQTL/T4.50000.cis_qtl_pairs.chr{autosome}.txt.gz", autosome = autosome)

    output:
        T0_nominal_merged = "results/eQTL/T0.50000.cis_qtl_pairs.ALL.txt.gz",
        T1_nominal_merged = "results/eQTL/T1.50000.cis_qtl_pairs.ALL.txt.gz",
        T2_nominal_merged = "results/eQTL/T2.50000.cis_qtl_pairs.ALL.txt.gz",
        T3_nominal_merged = "results/eQTL/T3.50000.cis_qtl_pairs.ALL.txt.gz",
        T4_nominal_merged = "results/eQTL/T4.50000.cis_qtl_pairs.ALL.txt.gz"

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
        T0_nominal_merged = "results/eQTL/T0.50000.cis_qtl_pairs.ALL.txt.gz",
        T1_nominal_merged = "results/eQTL/T1.50000.cis_qtl_pairs.ALL.txt.gz",
        T2_nominal_merged = "results/eQTL/T2.50000.cis_qtl_pairs.ALL.txt.gz",
        T3_nominal_merged = "results/eQTL/T3.50000.cis_qtl_pairs.ALL.txt.gz",
        T4_nominal_merged = "results/eQTL/T4.50000.cis_qtl_pairs.ALL.txt.gz"

    output:
        T0_unzipped_merged = "results/eQTL/T0.50000.cis_qtl_pairs.ALL.txt",
        T1_unzipped_merged = "results/eQTL/T1.50000.cis_qtl_pairs.ALL.txt",
        T2_unzipped_merged = "results/eQTL/T2.50000.cis_qtl_pairs.ALL.txt",
        T3_unzipped_merged = "results/eQTL/T3.50000.cis_qtl_pairs.ALL.txt",
        T4_unzipped_merged = "results/eQTL/T4.50000.cis_qtl_pairs.ALL.txt"

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
        permute = expand("results/eQTL/{timepoint}.50000.cis_qtl_fdr0.1.txt", timepoint = config["TIME"]),
        nominal = expand("results/eQTL/{timepoint}.50000.cis_qtl_pairs.ALL.txt", timepoint = config["TIME"]),
        script = "scripts/cis-eQTL_numbers.R"
    output:
        threshold_violin = "results/eQTL/results/violin_plot_pval_thresholds.pdf",
        eGene_tested = "results/eQTL/results/upset_cis_egenes_tested.pdf",
        upset_plot = "results/eQTL/results/upset_cis_egenes_fdr_0.1.pdf",
        stats = "results/eQTL/numbers_stats_cis.txt",
        TWAS_qtls = expand("results/eQTL/{timepoint}_cis_qtl_pairs_TWAS_input.txt", timepoint = config["TIME"])
    shell:
        '''
        Rscript {input.script} {input.permute} {input.nominal} {output.threshold_violin} {output.eGene_tested} {output.stats} {output.TWAS_qtls} {output.upset_plot}
        '''

