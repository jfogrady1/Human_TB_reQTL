autosomes = [str(i) for i in range(1, 23)]
rule all:
     input:
        #expand('work/RNA_seq/quantification/Normalised_counts_{timepoint}_qtltools.bed.gz', timepoint = config["TIME"]),
        #expand('data/covariate/{timepoint}_eqtlcovs_qtltools.txt', timepoint = config["TIME"]),
        #expand("work/RNA_seq/quantification/{timepoint}_residualised_expression.txt", timepoint = config["TIME"]),
        #DR_eQTLs = "results/reQTLs/DR_eQTLs_MASH_2_lfsr_0.05.txt",
        #"results/eQTL/results/INTERVAL_mashr_comparison.pdf"



rule edit_bed_file:
    input:
        expression_bed='work/RNA_seq/quantification/TMM_INT_counts_{timepoint}.bed.gz', # Phenotypes
        script = "scripts/bedfile_edit_qtltools.R" 
    output:
        bed = 'work/RNA_seq/quantification/Normalised_counts_{timepoint}_qtltools.bed.gz',
    
    params: 
        unzipped = 'work/RNA_seq/quantification/Normalised_counts_{timepoint}_qtltools.bed'
    shell:
        '''
        Rscript {input.script} {input.expression_bed} {params.unzipped} 
        bgzip {params.unzipped}
        tabix -p bed {output.bed}
        '''


rule edit_covariates:
    input:
        covariates_file='data/covariate/{timepoint}_eqtlcovs.txt', # Phenotypes
        script = "scripts/Covariate_edit_qtltools.R" 
    output:
        cov_edit = 'data/covariate/{timepoint}_eqtlcovs_qtltools.txt',
    
    shell:
        '''
        Rscript {input.script} {input.covariates_file} {output.cov_edit} 
        '''


rule extract_residuals:
    input:
        bed = multiext('work/RNA_seq/quantification/Normalised_counts_{timepoint}_qtltools.bed', ".gz", ".gz.tbi"), # Phenotypes
        covariates_file='data/covariate/{timepoint}_eqtlcovs_qtltools.txt'
    output:
        corrected = "work/RNA_seq/quantification/{timepoint}_residualised_expression.txt"
    singularity: "docker://jogrady/qtltools:1.3.1"
    shell:
        '''
        
        QTLtools correct --bed {input.bed} --cov {input.covariates_file} --out {output.corrected}
        '''

rule mashr:
    input:
        script = "scripts/MASHR.R",
        permute = expand("results/eQTL/{timepoint}.50000.cis_qtl_fdr0.1.txt", timepoint = config["TIME"]),
        symbols = "data/ref_genome/gencode.v43.annotation.gtf",
        corrected = expand("work/RNA_seq/quantification/{timepoint}_residualised_expression.txt", timepoint = config["TIME"]),
        vcf = "work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf.gz"

    output:
        significant_mash = "results/reQTLs/MashR_new_LFSR_values.txt",
        significant_mash_annotated = "results/reQTLs/MashR_new_LFSR_values_annotated.txt",
        DR_eQTLs = "results/reQTLs/DR_eQTLs_MASH_2_lfsr_0.05.txt",
        venn = "work/reQTLs/Venn_diagram_overlap_LFSR_0.05_MAG_2.pdf",
        t0_V_t1 = "work/reQTLs/T0_V_T1_reQTLs.txt",
        t0_V_t2 = "work/reQTLs/T0_V_T2_reQTLs.txt",
        t0_V_t3 = "work/reQTLs/T0_V_T3_reQTLs.txt",
        t0_V_t4 = "work/reQTLs/T0_V_T4_reQTLs.txt",
        common_genes = "work/reQTLs/Common_genes.txt",
        overlapping_reGenes = "work/reQTLs/Overlapping_response_eGenes_LFSR_0.05.txt",
        gprofiler_results = "work/reQTLs/Gprofiler_reGenes.txt",
        gprofiler_pdf = "results/reQTLs/Gprofiler_reGenes.pdf",
        mashr_plot_se = "results/reQTLs/MASHR_plot_effect_size.pdf",
        tss_response_eQTLs = "results/reQTLs/TSS_start_distance.pdf",
        TNFRSF10A = "results/reQTLs/TNFRS10A_response_eQTL.pdf",
        IFNGR2 = "results/reQTLs/IFNGR2_response_eQTL.pdf"


    params:
        time = expand("{t}", t = config["TIME"])

    shell:
        """
        Rscript {input.script} {params.time} {input.symbols} {output.significant_mash} {output.significant_mash_annotated} {output.DR_eQTLs} {output.venn} {output.t0_V_t1} \
        {output.t0_V_t2} {output.t0_V_t3} {output.t0_V_t4} {output.common_genes} {output.overlapping_reGenes} {output.gprofiler_results} {output.gprofiler_pdf} {output.mashr_plot_se} \
        {output.tss_response_eQTLs} {input.corrected} {input.vcf} {output.TNFRSF10A} {output.IFNGR2}
        """

rule eQTL_replication:
    input:
        script = "scripts/eQTL_replication.R",
        sig_results = expand("results/eQTL/{timepoint}.50000.cis_qtl_fdr0.1.txt", timepoint = config["TIME"]),
        nominal = expand("results/eQTL/{timepoint}.50000.cis_qtl_pairs.ALL.txt.gz", timepoint = config["TIME"]),
        INTERVAL_results = expand("data/gtex/analysis/INTERVAL_eQTL_summary_statistics/INTERVAL_eQTL_nominal_chr{chr}.tsv", chr = autosomes),
        INTERVAL_signif = "data/gtex/analysis/INTERVAL_eQTL_summary_statistics/INTERVAL_eQTL_phenotype_summary.tsv",
        mashr_results = "results/reQTLs/MashR_new_LFSR_values_annotated.txt",
    output:
        effect_comparison = "results/eQTL/results/Effect_size_slope.pdf",
        eQTL_correlation = "results/eQTL/results/INTERVAL_correaltion_comparison.pdf",
        mashr_correlation = "results/eQTL/results/INTERVAL_mashr_comparison.pdf"

    shell:
        """
        Rscript {input.script} {input.sig_results} \
        {output.effect_comparison} {input.nominal} \
        {input.INTERVAL_results} {input.INTERVAL_signif} {output.eQTL_correlation} {input.mashr_results} {output.mashr_correlation}
        """