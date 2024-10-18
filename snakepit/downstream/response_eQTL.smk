autosomes = [str(i) for i in range(1, 23)]
#rule all:
     #input:
        #expand('/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/Normalised_counts_{timepoint}_qtltools.bed.gz', timepoint = config["TIME"]),
        #expand('/home/workspace/jogrady/heQTL/data/covariate/{timepoint}_eqtlcovs_qtltools.txt', timepoint = config["TIME"]),
        #expand("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/{timepoint}_residualised_expression.txt", timepoint = config["TIME"]),
        #DR_eQTLs = "/home/workspace/jogrady/heQTL/results/reQTLs/DR_eQTLs_MASH_2_lfsr_0.05.txt",



rule edit_bed_file:
    input:
        expression_bed='/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/TMM_INT_counts_{timepoint}.bed.gz', # Phenotypes
        script = "/home/workspace/jogrady/heQTL/scripts/bedfile_edit_qtltools.R" 
    output:
        bed = '/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/Normalised_counts_{timepoint}_qtltools.bed.gz',
    
    params: 
        unzipped = '/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/Normalised_counts_{timepoint}_qtltools.bed'
    shell:
        '''
        Rscript {input.script} {input.expression_bed} {params.unzipped} 
        bgzip {params.unzipped}
        tabix -p bed {output.bed}
        '''


rule edit_covariates:
    input:
        covariates_file='/home/workspace/jogrady/heQTL/data/covariate/{timepoint}_eqtlcovs.txt', # Phenotypes
        script = "/home/workspace/jogrady/heQTL/scripts/Covariate_edit_qtltools.R" 
    output:
        cov_edit = '/home/workspace/jogrady/heQTL/data/covariate/{timepoint}_eqtlcovs_qtltools.txt',
    
    shell:
        '''
        Rscript {input.script} {input.covariates_file} {output.cov_edit} 
        '''


rule extract_residuals:
    input:
        bed = multiext('/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/Normalised_counts_{timepoint}_qtltools.bed', ".gz", ".gz.tbi"), # Phenotypes
        covariates_file='/home/workspace/jogrady/heQTL/data/covariate/{timepoint}_eqtlcovs_qtltools.txt'
    output:
        corrected = "/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/{timepoint}_residualised_expression.txt"
    singularity: "docker://jogrady/qtltools:1.3.1"
    shell:
        '''
        
        QTLtools correct --bed {input.bed} --cov {input.covariates_file} --out {output.corrected}
        '''

rule mashr:
    input:
        script = "/home/workspace/jogrady/heQTL/scripts/MASHR.R",
        permute = expand("/home/workspace/jogrady/heQTL/results/eQTL/{timepoint}.50000.cis_qtl_fdr0.1.txt", timepoint = config["TIME"]),
        symbols = "/home/workspace/jogrady/heQTL/data/ref_genome/gencode.v43.annotation.gtf",
        corrected = expand("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/{timepoint}_residualised_expression.txt", timepoint = config["TIME"]),
        vcf = "/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf.gz"

    output:
        significant_mash = "/home/workspace/jogrady/heQTL/results/reQTLs/MashR_new_LFSR_values.txt",
        significant_mash_annotated = "/home/workspace/jogrady/heQTL/results/reQTLs/MashR_new_LFSR_values_annotated.txt",
        DR_eQTLs = "/home/workspace/jogrady/heQTL/results/reQTLs/DR_eQTLs_MASH_2_lfsr_0.05.txt",
        venn = "/home/workspace/jogrady/heQTL/work/reQTLs/Venn_diagram_overlap_LFSR_0.05_MAG_2.pdf",
        t0_V_t1 = "/home/workspace/jogrady/heQTL/work/reQTLs/T0_V_T1_reQTLs.txt",
        t0_V_t2 = "/home/workspace/jogrady/heQTL/work/reQTLs/T0_V_T2_reQTLs.txt",
        t0_V_t3 = "/home/workspace/jogrady/heQTL/work/reQTLs/T0_V_T3_reQTLs.txt",
        t0_V_t4 = "/home/workspace/jogrady/heQTL/work/reQTLs/T0_V_T4_reQTLs.txt",
        common_genes = "/home/workspace/jogrady/heQTL/work/reQTLs/Common_genes.txt",
        overlapping_reGenes = "/home/workspace/jogrady/heQTL/work/reQTLs/Overlapping_response_eGenes_LFSR_0.05.txt",
        gprofiler_results = "/home/workspace/jogrady/heQTL/work/reQTLs/Gprofiler_reGenes.txt",
        gprofiler_pdf = "/home/workspace/jogrady/heQTL/results/reQTLs/Gprofiler_reGenes.pdf",
        mashr_plot_se = "/home/workspace/jogrady/heQTL/results/reQTLs/MASHR_plot_effect_size.pdf",
        tss_response_eQTLs = "/home/workspace/jogrady/heQTL/results/reQTLs/TSS_start_distance.pdf",
        TNFRSF10A = "/home/workspace/jogrady/heQTL/results/reQTLs/TNFRS10A_response_eQTL.pdf",
        IFNGR2 = "/home/workspace/jogrady/heQTL/results/reQTLs/IFNGR2_response_eQTL.pdf"


    params:
        time = expand("{t}", t = config["TIME"])

    shell:
        """
        Rscript {input.script} {params.time} {input.symbols} {output.significant_mash} {output.significant_mash_annotated} {output.DR_eQTLs} {output.venn} {output.t0_V_t1} \
        {output.t0_V_t2} {output.t0_V_t3} {output.t0_V_t4} {output.common_genes} {output.overlapping_reGenes} {output.gprofiler_results} {output.gprofiler_pdf} {output.mashr_plot_se} \
        {output.tss_response_eQTLs} {input.corrected} {input.vcf} {output.TNFRSF10A} {output.IFNGR2}
        """