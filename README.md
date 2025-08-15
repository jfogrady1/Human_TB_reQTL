# Genetic control of the transcriptional response to active tuberculosis disease and treatment

---

## Summary
The functional impact of genomic regulatory variation during active TB disease and anti-TB treatment (ATT) remains poorly understood. Here, using 240 paired whole-blood RNA-seq samples from n = 48 active TB patients who underwent ATT, we call and impute 1,506,948 genome-wide variants from these transcriptomes to map cis-expression quantitative trait loci (eQTL), response-eQTL, and cell type interaction eQTL (ieQTL) by computationally deconvolving the bulk RNA-seq data. We show that the whole-blood transcriptome is substantially perturbed during ATT by identifying 7,246 differentially expressed genes. We characterise 5,129 cis-eQTL and 587 response-eQTL and show that genes associated with reQTL are significantly enriched in pathways impacting the host response to mycobacterial challenge and xenobiotic metabolism. We highlight significant changes in cell type proportions due to ATT, notably for NK cells and classical monocytes. We characterise 921 ieQTL, including an active-TB classical monocyte-specific ieQTL for the gene NOD2. Our work sheds light on the immunogenetics of TB disease and treatment, while providing a framework for integrative genomics studies using only RNA-seq data.

---

## Overview of study

![Figure_01](https://github.com/user-attachments/assets/53a557b0-868f-4f99-8e53-5fc2a58c0fcb)

---

## Data Sources

* **RNA-seq Data:** The RNA-seq paired-end fastq files for n = 48 individuals (240 files in total) are accessible through the **Sequence Read Archive (SRA)** with the GEO accession number **GSE157657**.
* **scRNA-seq Data:** The PBMC single-cell RNA-seq data from n = 3 active TB patients and n = 2 LTBI patients are available on the **SRA** under accession number **SRP247583**. The specific SRR IDs are SRR11038989 through SRR11038993.
* **GIAB Datasets:** The Genome in a Bottle (GIAB) datasets can be accessed at this address: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release`.

### Reference Files and Software

* **Reference Genome:** The reference genome used for RNA-seq read alignment is **GRCh38_full_analysis_set_plus_decoy_hla**.
* **DeepVariant:** The **DeepVariant** variant calling model can be found on its official GitHub repository: `https://github.com/google/deepvariant`.

### Additional Data

Filtered and imputed variants for all n = 48 samples, as well as the raw cis-eQTL, reQTL, and ieQTL results, will be made publicly available upon publication of this study.

---
## Citation

O’Grady, J.F., Leonard, A.S., Li, H., Fang, L., Pausch, H., Gormley, I.C., Gordon, S.V., MacHugh, D.E. Genetic control of the transcriptional response to active tuberculosis disease and treatment [Preprint]. (2025) bioRxiv, 2025.04.11.648383. https://doi.org/10.1101/2025.04.11.648383. 
