# Genetic control of the transcriptional response to active tuberculosis disease and treatment

# Summary
Understanding the functional impact of variants is critical for discerning the host response during tuberculosis (TB) disease, and anti-TB treatment. Hitherto, there have been no genome-wide in vivo response expression quantitative trait loci (re-QTLs) in the context of active-TB and anti-TB treatment. Here, we exploit longitudinal whole-blood RNA-seq data from n = 48 patients with active-TB who underwent anti-TB treatment and call variants from these transcriptomics data. We highlight substantial differences in the transcriptome profile following anti-TB treatment with thousands of genes differentially expressed. Associating our imputed variants to the expression of nearby genes, we characterise thousands of cis-eQTLs and hundreds of re-QTLs. We further show significant changes in cell-type proportions during anti-TB treatment through in silico deconvolution and highlight the cell-type specific nature of cis-eQTLs. Our work highlights the immunobiology of active-TB and anti-TB treatment while also providing a framework for other studies that are limited to RNA-seq data.


## Overview of study

![Figure_01](https://github.com/user-attachments/assets/53a557b0-868f-4f99-8e53-5fc2a58c0fcb)

---

### Data Sources

* **RNA-seq Data:** The RNA-seq paired-end fastq files for n = 48 individuals (240 files in total) are accessible through the **Sequence Read Archive (SRA)** with the GEO accession number **GSE157657**.
* **scRNA-seq Data:** The PBMC single-cell RNA-seq data from n = 3 active TB patients and n = 2 LTBI patients are available on the **SRA** under accession number **SRP247583**. The specific SRR IDs are SRR11038989 through SRR11038993.
* **GIAB Datasets:** The Genome in a Bottle (GIAB) datasets can be accessed at this address: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release`.

### Reference Files and Software

* **Reference Genome:** The reference genome used for RNA-seq read alignment is **GRCh38_full_analysis_set_plus_decoy_hla**.
* **DeepVariant:** The **DeepVariant** variant calling model can be found on its official GitHub repository: `https://github.com/google/deepvariant`.

---

### Additional Data

Filtered and imputed variants for all n = 48 samples, as well as the raw cis-eQTL, reQTL, and ieQTL results, will be made publicly available upon publication of this study.
