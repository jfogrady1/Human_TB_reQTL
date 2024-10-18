#rule all:
     #input:
        #expand('/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/TMM_INT_counts_{timepoint}.bed.gz', timepoint = config["TIME"]),
        #expand('/home/workspace/jogrady/heQTL/work/ieQTL/{t}_NK_cell_interaction_input_from_cibersort.txt', t = config["TIME"], celltype = config["cell"]),
        #expand('/home/workspace/jogrady/heQTL/work/ieQTL/{t}_NK_cell.cis_qtl_top_assoc.txt.gz', t = config["TIME"], celltype = config["cell"]),
        #expand("/home/workspace/jogrady/heQTL/work/ieQTL/{t}_Memory_CD8+_T_cell.cis_qtl_top_assoc.txt.gz", t = config["TIME"]),
        #expand("/home/workspace/jogrady/heQTL/work/ieQTL/{t}_Naïve_B_cell.cis_qtl_top_assoc.txt.gz", t = config["TIME"]),
        #expand("/home/workspace/jogrady/heQTL/work/ieQTL/{t}_Classical_Monocyte.cis_qtl_top_assoc.txt.gz", t = config["TIME"]),
        #expand("/home/workspace/jogrady/heQTL/work/ieQTL/{t}_Non-classical_Monocyte.cis_qtl_top_assoc.txt.gz", t = config["TIME"]),
        #expand("/home/workspace/jogrady/heQTL/work/ieQTL/{t}_{celltype}_interaction_input_from_cibersort.txt", t = config["TIME"], celltype = config["cell"]),
        #expand("/home/workspace/jogrady/heQTL/work/ieQTL/{t}_cell_interaction_input_from_cibersort.txt", t = config["TIME"]),
        #'/home/workspace/jogrady/heQTL/results/ieQTLs/ieQTL_MASHR.RData'


rule preprocess_ieQTL:
    input:
        ciber = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/CIBERSORTx_Results_Real_deconvolution.txt",
        script = "/home/workspace/jogrady/heQTL/scripts/ieQTL_preprocessing_cibersort_files.R"
    output:
        cell_common = expand("/home/workspace/jogrady/heQTL/work/ieQTL/{t}_{celltype}_interaction_input_from_cibersort.txt", t = config["TIME"], celltype = config["cell"]),
        transformed = expand("/home/workspace/jogrady/heQTL/work/ieQTL/{t}_cell_interaction_input_from_cibersort.txt", t = config["TIME"])
    
    params:
        cell = expand("{celltype}", celltype = config["cell"]),
    resources:
        mem_mb = 6000,
        threads = 40
    shell:
        '''
        Rscript {input.script} {input.ciber} {params.cell} {output.transformed}
        '''      


rule interaction_NK:
    input:
        expression_bed='/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/TMM_INT_counts_{timepoint}.bed.gz', # Phenotypes
        covariates_file='/home/workspace/jogrady/heQTL/data/covariate/{timepoint}_eqtlcovs.txt', # Covariates,
        interaction_file = '/home/workspace/jogrady/heQTL/work/ieQTL/{timepoint}_NK_cell_interaction_input_from_cibersort.txt'
    output:
        parquet_files= "/home/workspace/jogrady/heQTL/work/ieQTL/{timepoint}_NK_cell.cis_qtl_top_assoc.txt.gz",
    resources:
        mem_mb = 6000,
        threads = 40
    params:
       plink_prefix_path="/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", # Genotypes
       prefix = "/home/workspace/jogrady/heQTL/work/ieQTL/{timepoint}_NK_cell",
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
        expression_bed='/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/TMM_INT_counts_{timepoint}.bed.gz', # Phenotypes
        covariates_file='/home/workspace/jogrady/heQTL/data/covariate/{timepoint}_eqtlcovs.txt', # Covariates,
        interaction_file = '/home/workspace/jogrady/heQTL/work/ieQTL/{timepoint}_Naïve_B_cell_interaction_input_from_cibersort.txt'
    output:
        parquet_files= "/home/workspace/jogrady/heQTL/work/ieQTL/{timepoint}_Naïve_B_cell.cis_qtl_top_assoc.txt.gz",
    resources:
        mem_mb = 6000,
        threads = 40
    params:
       plink_prefix_path="/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", # Genotypes
       prefix = "/home/workspace/jogrady/heQTL/work/ieQTL/{timepoint}_Naïve_B_cell",
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
        expression_bed='/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/TMM_INT_counts_{timepoint}.bed.gz', # Phenotypes
        covariates_file='/home/workspace/jogrady/heQTL/data/covariate/{timepoint}_eqtlcovs.txt', # Covariates,
        interaction_file = '/home/workspace/jogrady/heQTL/work/ieQTL/{timepoint}_Classical_Monocyte_interaction_input_from_cibersort.txt'
    output:
        parquet_files= "/home/workspace/jogrady/heQTL/work/ieQTL/{timepoint}_Classical_Monocyte.cis_qtl_top_assoc.txt.gz",
    resources:
        mem_mb = 6000,
        threads = 40
    params:
       plink_prefix_path="/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", # Genotypes
       prefix = "/home/workspace/jogrady/heQTL/work/ieQTL/{timepoint}_Classical_Monocyte",
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
        expression_bed='/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/TMM_INT_counts_{timepoint}.bed.gz', # Phenotypes
        covariates_file='/home/workspace/jogrady/heQTL/data/covariate/{timepoint}_eqtlcovs.txt', # Covariates,
        interaction_file = '/home/workspace/jogrady/heQTL/work/ieQTL/{timepoint}_Non-classical_Monocyte_interaction_input_from_cibersort.txt'
    output:
        parquet_files= "/home/workspace/jogrady/heQTL/work/ieQTL/{timepoint}_Non-classical_Monocyte.cis_qtl_top_assoc.txt.gz",
    resources:
        mem_mb = 6000,
        threads = 40
    params:
       plink_prefix_path="/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", # Genotypes
       prefix = "/home/workspace/jogrady/heQTL/work/ieQTL/{timepoint}_Non-classical_Monocyte",
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
        expression_bed='/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/TMM_INT_counts_{timepoint}.bed.gz', # Phenotypes
        covariates_file='/home/workspace/jogrady/heQTL/data/covariate/{timepoint}_eqtlcovs.txt', # Covariates,
        interaction_file = '/home/workspace/jogrady/heQTL/work/ieQTL/{timepoint}_Memory_CD8+_T_cell_interaction_input_from_cibersort.txt'
    output:
        parquet_files= "/home/workspace/jogrady/heQTL/work/ieQTL/{timepoint}_Memory_CD8+_T_cell.cis_qtl_top_assoc.txt.gz",
    resources:
        mem_mb = 6000,
        threads = 40
    params:
       plink_prefix_path="/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", # Genotypes
       prefix = "/home/workspace/jogrady/heQTL/work/ieQTL/{timepoint}_Memory_CD8+_T_cell",
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

rule ieQTL_mashr:
    input:
        script = "/home/workspace/jogrady/heQTL/scripts/ieQTL_mashR.R",
        annotation = "/home/workspace/jogrady/heQTL/data/ref_genome/gencode.v43.annotation.gtf"
    output: 
        circ_plot = "/home/workspace/jogrady/heQTL/work/ieQTL/MASHR_LFSR_0.01_5_Cell_types.pdf",
        all_signif_ieQTLs_long = "/home/workspace/jogrady/heQTL/results/ieQTLs/LFSR_0.01_ieQTLs.txt",
        private_ieQTL = "/home/workspace/jogrady/heQTL/results/ieQTLs/mTB_private_ieQTLs.txt",
        architecture = "/home/workspace/jogrady/heQTL/results/ieQTLs/Genomic_architecture_ieQTLs.pdf",
        rds = '/home/workspace/jogrady/heQTL/results/ieQTLs/ieQTL_MASHR.RData'

    params:
        cell_type = expand("{celltype}", celltype = config["cell"])

    shell:
        """
        Rscript {input.script} {params.cell_type} {output.circ_plot} {input.annotation} {output.all_signif_ieQTLs_long} {output.private_ieQTL} {output.architecture} {output.rds}
        """