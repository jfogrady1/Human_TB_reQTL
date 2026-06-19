# Genetic control of the transcriptional response to active tuberculosis disease and treatment

**Background:** The functional impact of genomic variation during active tuberculosis disease and anti-tuberculosis treatment, respectively, remains poorly understood owing to a paucity of in vivo context-specific genome-wide molecular quantitative trait loci studies for these phenotypic traits. Characterising the context-specific functional impact of sequence polymorphisms represents a necessary endeavour to provide new insight into how regulatory variation modulates complex traits, such as susceptibility to Mycobacterium tuberculosis infection.
**Methods:** Using 240 paired whole-blood RNA-seq samples from n = 48 active tuberculosis patients who subsequently underwent anti-TB treatment, we called and imputed genome-wide variants from these transcriptomes and mapped cis-expression quantitative trait loci (eQTL), response-eQTL (reQTL), and cell type interaction eQTL (ieQTL) after computationally deconvolving the bulk RNA-seq data using peripheral blood mononuclear single-cell RNA-seq data from n = 5 TB-infected patients as a reference. 
**Results:** Here, we characterise 1,506,948 high-quality imputed genome-wide variants from the transcriptomics data.  We identify a total of 5,356 cis-eQTL and 790 reQTL and show that genes associated with reQTL are significantly enriched in pathways that impact the host response to mycobacterial challenge and xenobiotic metabolism. Additionally, we highlight significant changes in the proportions of computationally deconvolved cell types during anti-tuberculosis treatment, notably for natural killer cells and classical monocytes. Leveraging these deconvolved cell type proportions, we characterise a total of 1,098 ieQTL, including an active-TB classical monocyte-specific ieQTL for the gene ALOX5.
**Conclusions:** Our work sheds light on the immunogenetics of tuberculosis disease and treatment and provides a framework for integrative genomics studies using only RNA-seq data.

## Overview of study

![Figure_01](https://github.com/jfogrady1/Human_TB_reQTL/blob/main/Figure_01.png)

## Data Sources

* **RNA-seq Data:** The RNA-seq paired-end fastq files for n = 48 individuals (240 files in total) are accessible through the **Sequence Read Archive (SRA)** with the GEO accession number **GSE157657**.
* **scRNA-seq Data:** The PBMC single-cell RNA-seq data from n = 3 active TB patients and n = 2 LTBI patients are available on the **SRA** under accession number **SRP247583**. The specific SRR IDs are SRR11038989 through SRR11038993.
* **GIAB Datasets:** The Genome in a Bottle (GIAB) datasets can be accessed at this address: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release`.

### Reference Files and Software

* **Reference Genome:** The reference genome used for RNA-seq read alignment is **GRCh38_full_analysis_set_plus_decoy_hla**.
* **DeepVariant:** The **DeepVariant** variant calling model can be found on its official GitHub repository: `https://github.com/google/deepvariant`.

### Additional Data

Filtered and imputed variants for all n = 48 samples, as well as the raw cis-eQTL, reQTL, and ieQTL results are available on Zenodo (https://zenodo.org/records/15103705)

## Citation

O’Grady, J.F., Leonard, A.S., Li, H., Fang, L., Pausch, H., Gormley, I.C., Gordon, S.V., MacHugh, D.E. Genetic control of the transcriptional response to active tuberculosis disease and treatment [Preprint]. (2025) bioRxiv, 2025.04.11.648383. https://doi.org/10.1101/2025.04.11.648383. 
