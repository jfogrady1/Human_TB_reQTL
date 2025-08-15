#rule all:
  #   input:
        #expand('work/RNA_seq/quantification/TMM_INT_counts_{timepoint}.bed.gz', timepoint = config["TIME"]),
        #expand('work/ieQTL/{t}_NK_cell_interaction_input_from_cibersort.txt', t = config["TIME"], celltype = config["cell"]),
        #expand('work/ieQTL/{t}_NK_cell.cis_qtl_top_assoc.txt.gz', t = config["TIME"], celltype = config["cell"]),
        #expand("work/ieQTL/{t}_Memory_CD8+_T_cell.cis_qtl_top_assoc.txt.gz", t = config["TIME"]),
        #expand("work/ieQTL/{t}_Naïve_B_cell.cis_qtl_top_assoc.txt.gz", t = config["TIME"]),
        #expand('data/covariate/{timepoint}_ieqtlcovs_{celltype}.txt', timepoint = config["TIME"], celltype = config["cell"]),
        #expand("work/ieQTL/{t}_Classical_Monocyte.cis_qtl_top_assoc.txt.gz", t = config["TIME"]),
        #expand("work/ieQTL/{t}_Non-classical_Monocyte.cis_qtl_top_assoc.txt.gz", t = config["TIME"]),
        #expand("work/ieQTL/{t}_{celltype}_interaction_input_from_cibersort.txt", t = config["TIME"], celltype = config["cell"]),
        #expand("work/ieQTL/{t}_cell_interaction_input_from_cibersort.txt", t = config["TIME"]),
        #expand("work/RNA_seq/quantification/{timepoint}_{celltype}_residualised_expression.txt", timepoint = config["TIME"], celltype = config["cell"]),
        #'results/ieQTLs/ieQTL_MASHR.RData'




rule preprocess_ieQTL:
    input:
        ciber = "work/scRNA_seq/Benchmark/CIBERSORTx_Results_Real_deconvolution.txt",
        script = "scripts/ieQTL_preprocessing_cibersort_files.R"
    output:
        cell_common = expand("work/ieQTL/{t}_{celltype}_interaction_input_from_cibersort.txt", t = config["TIME"], celltype = config["cell"]),
        transformed = expand("work/ieQTL/{t}_cell_interaction_input_from_cibersort.txt", t = config["TIME"])
    
    params:
        cell = expand("{celltype}", celltype = config["cell"]),
    resources:
        mem_mb = 6000,
        threads = 40
    shell:
        '''
        Rscript {input.script} {input.ciber} {params.cell} {output.transformed}
        '''      

rule add_celltype_to_covariates:
    input:
        script = "scripts/ieQTL_covariates_preprocessing.R"
    output:
        covariate_files = expand('data/covariate/{timepoint}_ieqtlcovs_{celltype}.txt', timepoint = config["TIME"], celltype = config["cell"])
    shell:
        """
        Rscript {input.script}
        """

rule interaction_NK:
    input:
        expression_bed='work/RNA_seq/quantification/TMM_INT_counts_{timepoint}.bed.gz', # Phenotypes
        covariates_file='data/covariate/{timepoint}_ieqtlcovs_NK_cell.txt', # Covariates,
        interaction_file = 'work/ieQTL/{timepoint}_NK_cell_interaction_input_from_cibersort.txt'
    output:
        parquet_files= "work/ieQTL/{timepoint}_NK_cell.cis_qtl_top_assoc.txt.gz",
    resources:
        mem_mb = 6000,
        threads = 40
    params:
       plink_prefix_path="work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", # Genotypes
       prefix = "work/ieQTL/{timepoint}_NK_cell",
       cis_window_size = "50000"
    shell:
        '''
        python3 -m tensorqtl {params.plink_prefix_path} {input.expression_bed} {params.prefix} \
            --covariates {input.covariates_file} \
            --window {params.cis_window_size} \
            --interaction {input.interaction_file} \
            --maf_threshold_interaction 0.05 \
            --mode cis_nominal
        '''


rule interaction_Naive_B:
    input:
        expression_bed='work/RNA_seq/quantification/TMM_INT_counts_{timepoint}.bed.gz', # Phenotypes
        covariates_file='data/covariate/{timepoint}_ieqtlcovs_Naïve_B_cell.txt', # Covariates,
        interaction_file = 'work/ieQTL/{timepoint}_Naïve_B_cell_interaction_input_from_cibersort.txt'
    output:
        parquet_files= "work/ieQTL/{timepoint}_Naïve_B_cell.cis_qtl_top_assoc.txt.gz",
    resources:
        mem_mb = 6000,
        threads = 40
    params:
       plink_prefix_path="work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", # Genotypes
       prefix = "work/ieQTL/{timepoint}_Naïve_B_cell",
       cis_window_size = "50000"
    shell:
        '''
        python3 -m tensorqtl {params.plink_prefix_path} {input.expression_bed} {params.prefix} \
            --covariates {input.covariates_file} \
            --window {params.cis_window_size} \
            --interaction {input.interaction_file} \
            --maf_threshold_interaction 0.05 \
            --mode cis_nominal
        '''
rule interaction_Monocyte:
    input:
        expression_bed='work/RNA_seq/quantification/TMM_INT_counts_{timepoint}.bed.gz', # Phenotypes
        covariates_file='data/covariate/{timepoint}_ieqtlcovs_Classical_Monocyte.txt', # Covariates,
        interaction_file = 'work/ieQTL/{timepoint}_Classical_Monocyte_interaction_input_from_cibersort.txt'
    output:
        parquet_files= "work/ieQTL/{timepoint}_Classical_Monocyte.cis_qtl_top_assoc.txt.gz",
    resources:
        mem_mb = 6000,
        threads = 40
    params:
       plink_prefix_path="work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", # Genotypes
       prefix = "work/ieQTL/{timepoint}_Classical_Monocyte",
       cis_window_size = "50000"
    shell:
        '''
        python3 -m tensorqtl {params.plink_prefix_path} {input.expression_bed} {params.prefix} \
            --covariates {input.covariates_file} \
            --window {params.cis_window_size} \
            --interaction {input.interaction_file} \
            --maf_threshold_interaction 0.05 \
            --mode cis_nominal
        '''
rule interaction_non_classical_monocyte:
    input:
        expression_bed='work/RNA_seq/quantification/TMM_INT_counts_{timepoint}.bed.gz', # Phenotypes
        covariates_file='data/covariate/{timepoint}_ieqtlcovs_Non-classical_Monocyte.txt', # Covariates,
        interaction_file = 'work/ieQTL/{timepoint}_Non-classical_Monocyte_interaction_input_from_cibersort.txt'
    output:
        parquet_files= "work/ieQTL/{timepoint}_Non-classical_Monocyte.cis_qtl_top_assoc.txt.gz",
    resources:
        mem_mb = 6000,
        threads = 40
    params:
       plink_prefix_path="work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", # Genotypes
       prefix = "work/ieQTL/{timepoint}_Non-classical_Monocyte",
       cis_window_size = "50000"
    shell:
        '''
        python3 -m tensorqtl {params.plink_prefix_path} {input.expression_bed} {params.prefix} \
            --covariates {input.covariates_file} \
            --window {params.cis_window_size} \
            --interaction {input.interaction_file} \
            --maf_threshold_interaction 0.05 \
            --mode cis_nominal
        '''

rule interaction_CD8_Tcell:
    input:
        expression_bed='work/RNA_seq/quantification/TMM_INT_counts_{timepoint}.bed.gz', # Phenotypes
        covariates_file='data/covariate/{timepoint}_ieqtlcovs_Memory_CD8+_T_cell.txt', # Covariates,
        interaction_file = 'work/ieQTL/{timepoint}_Memory_CD8+_T_cell_interaction_input_from_cibersort.txt'
    output:
        parquet_files= "work/ieQTL/{timepoint}_Memory_CD8+_T_cell.cis_qtl_top_assoc.txt.gz",
    resources:
        mem_mb = 6000,
        threads = 40
    params:
       plink_prefix_path="work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", # Genotypes
       prefix = "work/ieQTL/{timepoint}_Memory_CD8+_T_cell",
       cis_window_size = "50000"
    shell:
        '''
        python3 -m tensorqtl {params.plink_prefix_path} {input.expression_bed} {params.prefix} \
            --covariates {input.covariates_file} \
            --window {params.cis_window_size} \
            --interaction {input.interaction_file} \
            --maf_threshold_interaction 0.05 \
            --mode cis_nominal
        '''




wildcard_constraints:
    timepoint="T[0-9]+",  # Matches T0, T1, T2, etc.
    celltype="[A-Za-z0-9+_ïë-]+"
       
rule edit_ieQTL_covariates:
    input:
        covariates_file='data/covariate/{timepoint}_ieqtlcovs_{celltype}.txt',
        script="scripts/Covariate_edit_qtltools.R"
    output:
        cov_edit='data/covariate/{timepoint}_ieqtlcovs_{celltype}_qtltools.txt'
    shell:
        '''
        Rscript {input.script} {input.covariates_file} {output.cov_edit} 
        '''
        
rule extract_residuals:
    input:
        bed = multiext('work/RNA_seq/quantification/Normalised_counts_{timepoint}_qtltools.bed', ".gz", ".gz.tbi"), # Phenotypes
        covariates_file='data/covariate/{timepoint}_ieqtlcovs_{celltype}_qtltools.txt'
    output:
        corrected = "work/RNA_seq/quantification/{timepoint}_{celltype}_residualised_expression.txt"
    singularity: "docker://jogrady/qtltools:1.3.1"
    shell:
        '''
        
        QTLtools correct --bed {input.bed} --cov {input.covariates_file} --out {output.corrected}
        '''

rule ieQTL_mashr:
    input:
        script = "scripts/ieQTL_mashR.R",
        annotation = "data/ref_genome/gencode.v43.annotation.gtf",
        ieQTL= expand("work/ieQTL/{timepoint}_{celltype}.cis_qtl_top_assoc.txt.gz", timepoint = config["TIME"], celltype = config["cell"]), # note not used at all here

    output: 
        circ_plot = "work/ieQTL/MASHR_LFSR_0.01_5_Cell_types.pdf",
        all_signif_ieQTLs_long = "results/ieQTLs/LFSR_0.01_ieQTLs.txt",
        private_ieQTL = "results/ieQTLs/mTB_private_ieQTLs.txt",
        private_treat = "results/ieQTLs/treat_private_ieQTLs.txt",
        architecture = "results/ieQTLs/Genomic_architecture_ieQTLs.pdf",
        rds = 'results/ieQTLs/ieQTL_MASHR.RData'

    params:
        cell_type = expand("{celltype}", celltype = config["cell"])

    shell:
        """
        Rscript {input.script} {params.cell_type} {output.circ_plot} {input.annotation} {output.all_signif_ieQTLs_long} {output.private_ieQTL} {output.private_treat} {output.architecture} {output.rds}
        """