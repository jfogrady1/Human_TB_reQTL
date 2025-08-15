#rule all:
 #  input:
        #'work/DNA_seq/imputation/final/chr14_karyogram_raw.pdf',
        #'work/DNA_seq/imputation/final/PCA_pruned_imputed_variants.pdf',
        #'work/DNA_seq/imputation/final/Chr14_CDS.bed',
        #'work/DNA_seq/imputation/final/Chr14.vcf',
        #'work/DNA_seq/imputation/final/Typed_variants_chr14.bed',
        #'work/DNA_seq/imputation/final/Imputed_variants_chr14.bed',
        #'work/DNA_seq/imputation/final/Imputed_variants_R0.3_chr14.bed',
        #expand('work/DNA_seq/imputation/final/Overlap_Chr14_CDS_R0.{R}.txt', R = [6,3]),
        #'work/DNA_seq/imputation/final/SNP_density_all_chr.pdf',
        #'work/DNA_seq/imputation/final/chr14_karyogram_raw.pdf'

rule LD_prune:
    input:
        imputed = 'work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf.gz'
    output:
        pruned = 'work/DNA_seq/imputation/final/DV_variants.prune.in',
        eigenvec =   'work/DNA_seq/imputation/final/DV_pruned.eigenvec',
        eigenval = 'work/DNA_seq/imputation/final/DV_pruned.eigenval'
    params:
        prefix = "work/DNA_seq/imputation/final/DV_variants",
        prefix_pruned = "work/DNA_seq/imputation/final/DV_pruned"       
    shell:
        '''
        plink --vcf {input.imputed} --keep-allele-order --autosome --indep-pairwise 50 5 0.01 --out {params.prefix}

        plink --vcf {input.imputed} --keep-allele-order --autosome --extract {params.prefix}.prune.in --make-bed --pca --out {params.prefix_pruned}
        '''
rule PCA:
    input:
        eigenvec =   'work/DNA_seq/imputation/final/DV_pruned.eigenvec',
        eigenval = 'work/DNA_seq/imputation/final/DV_pruned.eigenval',
        ancestry = 'data/covariate/Ancestry.txt',
        script = 'scripts/genotype_PCA.R'
    output:
        pca = 'work/DNA_seq/imputation/final/PCA_pruned_imputed_variants.pdf'
    
    shell:
        '''
        Rscript {input.script} {input.eigenvec} {input.eigenval} {input.ancestry} {output.pca}
        '''

rule Chr14_extract:
    input:
         raw = 'work/DNA_seq/imputation/final/DV_raw_IMPUTED.vcf.gz'
    output:
        chr14 = 'work/DNA_seq/imputation/final/chr14_dr2.txt'
    shell:
        '''
         software/bcftools-1.15.1/bcftools view {input.raw} --regions chr14 | software/bcftools-1.15.1/bcftools query -f '%CHROM\t%POS\t%INFO\n' > {output.chr14}
        '''

rule Chr14_plotting:
    input:
        counts = 'work/RNA_seq/quantification/count_matrix_ordered_ALL.txt',
        annotation = 'data/ref_genome/gencode.v43.annotation.gtf',
        SNP_data = 'work/DNA_seq/imputation/final/chr14_dr2.txt',
        script = 'scripts/Chr14_karyoplot.R'    
    output:
        karyogram = 'work/DNA_seq/imputation/final/chr14_karyogram_raw.pdf',        
    shell:
        '''
        Rscript {input.script} {input.counts} {input.annotation} {input.SNP_data} {output.karyogram}
        '''



rule Chr14_CDS_extract:
    input:
        annotation = 'data/GIAB/gencode.v43.basic.annotation.gff3.gz',
    output:
        cds = 'work/DNA_seq/imputation/final/Chr14_CDS.bed'
    shell:
        '''
        gzip -dc {input.annotation} | awk -v OFS='\t' '$1 == "chr14" && $3 == "gene"  && $4 < $5 {{ print $1, $4, $5, "CDS" }}' | awk '!dup[$0]++' > {output.cds}
        '''

rule Chr14_vcf:
     input:
        vcf = 'work/DNA_seq/imputation/final/DV_raw_IMPUTED.vcf.gz',
        vcf60 = 'work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf.gz',
        vcf30 = 'work/DNA_seq/imputation/final/DV_IMPUTED_R0.3_filtered.vcf.gz'
     output:
        unzipped = 'work/DNA_seq/imputation/final/Chr14.vcf',
        unzipped60 = 'work/DNA_seq/imputation/final/Chr14_R0.6.vcf',
        unzipped30 = 'work/DNA_seq/imputation/final/Chr14_R0.3.vcf'
     shell:
        '''
        software/bcftools-1.15.1/bcftools view {input.vcf} --regions chr14 -Ov -o {output.unzipped}
        software/bcftools-1.15.1/bcftools index -f {input.vcf60}
        software/bcftools-1.15.1/bcftools view {input.vcf60} --regions chr14 -Ov -o {output.unzipped60}
        software/bcftools-1.15.1/bcftools index -f {input.vcf30}
        software/bcftools-1.15.1/bcftools view {input.vcf30} --regions chr14 -Ov -o {output.unzipped30}
        '''

rule vcf2bed:
     input:
        vcf = 'work/DNA_seq/imputation/final/Chr14.vcf',
        vcf60 = 'work/DNA_seq/imputation/final/Chr14_R0.6.vcf',
        vcf30 = 'work/DNA_seq/imputation/final/Chr14_R0.3.vcf'

     output:
        typed = 'work/DNA_seq/imputation/final/Typed_variants_chr14.bed',
        imputed = 'work/DNA_seq/imputation/final/Imputed_variants_chr14.bed',
        imputed60 = 'work/DNA_seq/imputation/final/Imputed_variants_R0.6_chr14.bed',
        imputed30 = 'work/DNA_seq/imputation/final/Imputed_variants_R0.3_chr14.bed'
     
     shell:
        '''
        # load bedops
        vcf2bed < {input.vcf} | grep -v "IMP" | cut -f 1-3 > {output.typed}

        vcf2bed < {input.vcf} | grep "IMP" | cut -f 1-3 > {output.imputed}

        vcf2bed < {input.vcf60} | grep "IMP" | cut -f 1-3 > {output.imputed60}

        vcf2bed < {input.vcf30} | grep "IMP" | cut -f 1-3 > {output.imputed30}
        '''

rule overlap:
    input:
        typed = 'work/DNA_seq/imputation/final/Typed_variants_chr14.bed',
        imputed = 'work/DNA_seq/imputation/final/Imputed_variants_R0.{R}_chr14.bed',
        cds = 'work/DNA_seq/imputation/final/Chr14_CDS.bed'
    output:  
        overlap = 'work/DNA_seq/imputation/final/Overlap_Chr14_CDS_R0.{R}.txt'
        
    shell:
        '''
        echo "Total TYPED" > {output.overlap}
        wc -l {input.typed} >> {output.overlap}
        echo "Total TYPED (In GENES)" >> {output.overlap}
        bedtools intersect -a {input.typed} -b {input.cds} -u | wc -l >> {output.overlap}

        echo "Total IMPUTED" >> {output.overlap}
        wc -l {input.imputed} >> {output.overlap}
        echo "Total IMPUTED (In GENES)" >> {output.overlap}
        bedtools intersect -a {input.imputed} -b {input.cds} -u | wc -l >> {output.overlap}
        '''


rule density_plot:
     input:
        vcf = 'work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf.gz',
        script = 'scripts/SNP_density_plotting.R'
     output:
        plot = 'work/DNA_seq/imputation/final/SNP_density_all_chr.pdf',
        zoom1 = 'work/DNA_seq/imputation/final/SNP_density_chr11_chr15.pdf',
        zoom2 = 'work/DNA_seq/imputation/final/SNP_density_chr13_zoom.pdf',
     shell:
        '''
        software/bcftools-1.15.1/bcftools query -f '%ID\t%CHROM\t%POS\t%POS\t%INFO\n' {input.vcf} | sed 's/chr//g' > work/DNA_seq/imputation/final/snp_density.txt

        Rscript {input.script} work/DNA_seq/imputation/final/snp_density.txt {output.plot} {output.zoom1} {output.zoom2} 
        '''
rule vep:
    input:
        vep = 'software/enembl-vep/vep',
        vcf = 'work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf.gz',

    output:
        vep = 'software/enembl-vep/VEP_human_analysis.txt'
    shell:
        '''
        {input.vep} -i  {input.vcf}  -o VEP_human_analysis.txt --offline -force_overwrite --canonical --symbol --tab --pick --fields Uploaded_variation,SYMBOL,Feature,Feture_type,Consequence,cDNA_position,CANONICAL
        '''
rule vep_plot:
    input:
        vcf = 'work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf.gz',
        script = 'scripts/VEP_analysis.R',
        vep_output = 'software/enembl-vep/VEP_human_analysis.txt'
    output:
        plot = "results/eQTL/results/Variant_location_called_V_Imputed.pdf"