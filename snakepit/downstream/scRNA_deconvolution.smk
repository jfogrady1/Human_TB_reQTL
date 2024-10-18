import pandas as pd
import os
import subprocess

#rule all:
    #input:
        #expand('/home/workspace/jogrady/heQTL/work/scRNA_seq/{singlecell_ids}', singlecell_ids=config["sc_ids"]),
        #'/home/workspace/jogrady/heQTL/work/scRNA_seq/TB.combined.metadata.txt',
        #'/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/ALL_benchmark_results.txt',
        #"/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/Raw_count_matrix_all_samples.txt",
        #"/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/CIBERSORTx_Music_correlation_res.txt"

rule cellranger_count:
    input:
        reference = "/home/workspace/jogrady/heQTL/data/scRNA_reference/refdata-gex-GRCh38-2020-A",
        fastq_directory = "/home/workspace/jogrady/heQTL/data/scRNA_reference/",
        cell_ranger = "/home/workspace/jogrady/heQTL/software/cellranger-7.2.0/cellranger"

    output:
        '/home/workspace/jogrady/heQTL/work/scRNA_seq/{singlecell_ids}'
    params:
        id_flag = "{singlecell_ids}"

    shell:
        """
         {input.cell_ranger} count --id {params.id_flag} --localcores 10 --expect-cells 10000 --no-bam --transcriptome {input.reference} --fastqs {input.fastq_directory} --sample {params.id_flag}
        """

#rule move_cell_ranger:
 #   input:
  #      counts = '/home/workspace/jogrady/heQTL/{singlecell_ids}'
#
 #   output:
  #      new_counts = directory('/home/workspace/jogrady/heQTL/work/scRNA_seq/{singlecell_ids}')
#
 #   shell:
  #      '''
   #     mv {input.counts} {output.new_counts}
    #    '''


rule scRNA_seq_preprocessing:
    input:
        sc_ids = expand('/home/workspace/jogrady/heQTL/work/scRNA_seq/{singlecell_ids}', singlecell_ids=config["sc_ids"]),
        SCINA_reference = '/home/workspace/jogrady/heQTL/data/scRNA_reference/pbmc_22_10x_cell_type_signature_gene_sets.gmt',
        script = "/home/workspace/jogrady/heQTL/scripts/scRNA_preprocessing.R"

    output:
        UMAP_no_annotation = "/home/workspace/jogrady/heQTL/work/scRNA_seq/UMAP_no_annotation_merged.pdf",
        TB_combined_seed = "/home/workspace/jogrady/heQTL/work/scRNA_seq/TB.combined.seed.set.RData",
        UMAP_annotated ="/home/workspace/jogrady/heQTL/work/scRNA_seq/UMAP_TB.integrated.annotated.pdf",
        TB_combined = "/home/workspace/jogrady/heQTL/work/scRNA_seq/TB.combined.final.rds",
        all_conserved_unique = '/home/workspace/jogrady/heQTL/work/scRNA_seq/Conserved_Marker_genes.txt',
        Active_id ='/home/workspace/jogrady/heQTL/work/scRNA_seq/Convert_UMI_Label.tsv',
        gene_counts_cell = '/home/workspace/jogrady/heQTL/work/scRNA_seq/Gene_Count_per_Cell.tsv',
        pseudobulk = '/home/workspace/jogrady/heQTL/work/scRNA_seq/TB.combined.pseudobulk.txt',
        metadata = '/home/workspace/jogrady/heQTL/work/scRNA_seq/TB.combined.metadata.txt'

    shell:
        """
        Rscript {input.script} {input.sc_ids} {input.SCINA_reference} {output.UMAP_no_annotation} {output.TB_combined_seed} \
        {output.UMAP_annotated} {output.TB_combined} {output.all_conserved_unique} {output.Active_id} {output.gene_counts_cell} \
        {output.pseudobulk} {output.metadata}
        """


rule scRNA_seq_deconvolution:
    input:
        annotation = "/home/workspace/jogrady/heQTL/data/ref_genome/gencode.v43.annotation.gtf",
        markers = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Conserved_Marker_genes.txt",
        sc_data = "/home/workspace/jogrady/heQTL/work/scRNA_seq/TB.combined.final.rds",
        Ciber_RAW_RAW = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/CIBERSORTx_RAW_RAW.txt",
        Ciber_RAW_TPM = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/CIBERSORTx_RAW_TPM.txt",
        Ciber_RAW_CPM = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/CIBERSORTx_RAW_CPM.txt",
        Ciber_TPM_RAW = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/CIBERSORTx_TPM_RAW.txt",
        Ciber_TPM_CPM = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/CIBERSORTx_TPM_CPM.txt",
        Ciber_TPM_TPM = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/CIBERSORTx_TPM_TPM.txt",
        script = "/home/workspace/jogrady/heQTL/scripts/scRNA_deconvolution_MuSiC.R",


    output:
        TPM_matrix = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/TPM_signature_matrix_Cibersort.txt", 
        RAW_matrix = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/RAW_signature_matrix_Cibersort.txt", 
        simbu_generation = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/Simbu_Simulation_Cell_Components.pdf",
        bulk_CPM = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/Pseudobulk_CPM.txt",
        bulk_TPM = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/Pseudobulk_TPM.txt",
        bulk_RAW = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/Pseudobulk_RAW.txt",
        RMSE_heatmap = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/Deconvolution_Benchmark_heatmap.pdf",
        Pearson_heatmap = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/Deconvolution_Benchmark_heatmap_pearson.pdf",
        Benchmark_no_ciber = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/Deconvolution_Pseudo_Benchmark_raw.pdf",
        Benchmark_ciber = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/Deconvolution_Pseudo_Benchmark_Cibersort_raw.pdf",
        Cibersort_CPM_individual = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/Benchmark_Cibersort_all_celltypes_RAW_CPM.pdf",
        All_benchmark_results = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/ALL_benchmark_results.txt"

    shell:
        """
        Rscript {input.script} {input.annotation} {input.markers} {input.sc_data} {output.TPM_matrix} {output.RAW_matrix} \
        {output.simbu_generation} {output.bulk_CPM} {output.bulk_TPM} {output.bulk_RAW} {input.Ciber_RAW_RAW} {input.Ciber_RAW_TPM} \
        {input.Ciber_RAW_CPM} {input.Ciber_TPM_RAW} {input.Ciber_TPM_TPM} {input.Ciber_TPM_CPM} {output.RMSE_heatmap} {output.Pearson_heatmap} \
        {output.Benchmark_no_ciber} {output.Benchmark_ciber} {output.Cibersort_CPM_individual} {output.All_benchmark_results}
        """


rule get_count_matrix_4_real_decon:
    input:
        counts = expand("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_{timepoint}.txt", timepoint = config["TIME"]),
        symbols ="/home/workspace/jogrady/heQTL/data/ref_genome/gencode.v43.annotation.gtf",
        script = "/home/workspace/jogrady/heQTL/scripts/Aggregate_counts_for_cibersort.R"
    output:
        merged = "/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/Raw_count_matrix_all_samples.txt"
    shell:
        """
        Rscript {input.script} {input.counts} {input.symbols} {output.merged}
        """


rule Cibersort_plotting:
    input:
        script = "/home/workspace/jogrady/heQTL/scripts/Cibersort_plotting.R",
        Ciber_decon_res = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/CIBERSORTx_Results_Real_deconvolution.txt",
        Raw_RNA_seq_counts = "/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/Raw_count_matrix_all_samples.txt",
        markers = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Conserved_Marker_genes.txt",
        tb_combined_file = "/home/workspace/jogrady/heQTL/work/scRNA_seq/TB.combined.final.rds"

    output:
        limma_results = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/Limma_results_real_deconvolution.txt",
        ciber_real_deconvolution = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/Cibersort_Deconvolved_real_data.pdf",
        music_ciber_compar = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/CIBERSORTx_Music_raw_comparison.pdf",
        correlation = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/CIBERSORTx_Music_correlation.pdf",
        correlation_table = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/CIBERSORTx_Music_correlation_res.txt"

    shell:
        """
        Rscript {input.script} {input.Ciber_decon_res} {output.limma_results} {output.ciber_real_deconvolution} {input.Raw_RNA_seq_counts} {input.markers} {input.tb_combined_file} {output.music_ciber_compar} \
        {output.correlation} {output.correlation_table}

        """


#rule ieQTL_preprocessing:
 #   # Generates INT transformed cell type proportions for all Cibersort values
    # Also tests for significant differences in the abundance of particular cell types
#
 #   input:
  #      Cibersort_output = "/home/workspace/jogrady/heQTL/work/scRNA_seq/CIBERSORTx_Results_Real_deconvolution.txt""
   # output:
    #    Cibersort_graph = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/Cibersort_Deconvolved_real_data.pdf"
#
#
 #   shell:
  #      """
#
#
 #       """
