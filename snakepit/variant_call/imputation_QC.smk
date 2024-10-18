#rule all:
 #  input:
        #'/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/chr14_karyogram_raw.pdf',
        #'/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/PCA_pruned_imputed_variants.pdf',
        #'/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Chr14_CDS.bed',
        #'/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Chr14.vcf',
        #'/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Typed_variants_chr14.bed',
        #'/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Imputed_variants_chr14.bed',
        #'/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Imputed_variants_R0.3_chr14.bed',
        #expand('/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Overlap_Chr14_CDS_R0.{R}.txt', R = [6,3]),
        #'/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/SNP_density_all_chr.pdf',
        #'/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/chr14_karyogram_raw.pdf'

rule LD_prune:
    input:
        imputed = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf.gz'
    output:
        pruned = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_variants.prune.in',
        eigenvec =   '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_pruned.eigenvec',
        eigenval = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_pruned.eigenval'
    params:
        prefix = "/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_variants",
        prefix_pruned = "/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_pruned"       
    shell:
        '''
        plink --vcf {input.imputed} --keep-allele-order --autosome --indep-pairwise 50 5 0.01 --out {params.prefix}

        plink --vcf {input.imputed} --keep-allele-order --autosome --extract {params.prefix}.prune.in --make-bed --pca --out {params.prefix_pruned}
        '''
rule PCA:
    input:
        eigenvec =   '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_pruned.eigenvec',
        eigenval = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_pruned.eigenval',
        ancestry = '/home/workspace/jogrady/heQTL/data/covariate/Ancestry.txt',
        script = '/home/workspace/jogrady/heQTL/scripts/genotype_PCA.R'
    output:
        pca = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/PCA_pruned_imputed_variants.pdf'
    
    shell:
        '''
        Rscript {input.script} {input.eigenvec} {input.eigenval} {input.ancestry} {output.pca}
        '''

rule Chr14_extract:
    input:
         raw = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_raw_IMPUTED.vcf.gz'
    output:
        chr14 = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/chr14_dr2.txt'
    shell:
        '''
         /home/workspace/jogrady/heQTL/software/bcftools-1.15.1/bcftools view {input.raw} --regions chr14 | /home/workspace/jogrady/heQTL/software/bcftools-1.15.1/bcftools query -f '%CHROM\t%POS\t%INFO\n' > {output.chr14}
        '''

rule Chr14_plotting:
    input:
        counts = '/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_ALL.txt',
        annotation = '/home/workspace/jogrady/heQTL/data/ref_genome/gencode.v43.annotation.gtf',
        SNP_data = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/chr14_dr2.txt',
        script = '/home/workspace/jogrady/heQTL/scripts/Chr14_karyoplot.R'    
    output:
        karyogram = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/chr14_karyogram_raw.pdf',        
    shell:
        '''
        Rscript {input.script} {input.counts} {input.annotation} {input.SNP_data} {output.karyogram}
        '''



rule Chr14_CDS_extract:
    input:
        annotation = '/home/workspace/jogrady/heQTL/data/GIAB/gencode.v43.basic.annotation.gff3.gz',
    output:
        cds = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Chr14_CDS.bed'
    shell:
        '''
        gzip -dc {input.annotation} | awk -v OFS='\t' '$1 == "chr14" && $3 == "gene"  && $4 < $5 {{ print $1, $4, $5, "CDS" }}' | awk '!dup[$0]++' > {output.cds}
        '''

rule Chr14_vcf:
     input:
        vcf = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_raw_IMPUTED.vcf.gz',
        vcf60 = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf.gz',
        vcf30 = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_IMPUTED_R0.3_filtered.vcf.gz'
     output:
        unzipped = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Chr14.vcf',
        unzipped60 = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Chr14_R0.6.vcf',
        unzipped30 = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Chr14_R0.3.vcf'
     shell:
        '''
        /home/workspace/jogrady/heQTL/software/bcftools-1.15.1/bcftools view {input.vcf} --regions chr14 -Ov -o {output.unzipped}
        /home/workspace/jogrady/heQTL/software/bcftools-1.15.1/bcftools index -f {input.vcf60}
        /home/workspace/jogrady/heQTL/software/bcftools-1.15.1/bcftools view {input.vcf60} --regions chr14 -Ov -o {output.unzipped60}
        /home/workspace/jogrady/heQTL/software/bcftools-1.15.1/bcftools index -f {input.vcf30}
        /home/workspace/jogrady/heQTL/software/bcftools-1.15.1/bcftools view {input.vcf30} --regions chr14 -Ov -o {output.unzipped30}
        '''

rule vcf2bed:
     input:
        vcf = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Chr14.vcf',
        vcf60 = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Chr14_R0.6.vcf',
        vcf30 = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Chr14_R0.3.vcf'

     output:
        typed = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Typed_variants_chr14.bed',
        imputed = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Imputed_variants_chr14.bed',
        imputed60 = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Imputed_variants_R0.6_chr14.bed',
        imputed30 = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Imputed_variants_R0.3_chr14.bed'
     
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
        typed = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Typed_variants_chr14.bed',
        imputed = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Imputed_variants_R0.{R}_chr14.bed',
        cds = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Chr14_CDS.bed'
    output:  
        overlap = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/Overlap_Chr14_CDS_R0.{R}.txt'
        
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
        vcf = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf.gz',
        script = '/home/workspace/jogrady/heQTL/scripts/SNP_density_plotting.R'
     output:
        plot = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/SNP_density_all_chr.pdf',
        zoom1 = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/SNP_density_chr11_chr15.pdf',
        zoom2 = '/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/SNP_density_chr13_zoom.pdf',
     shell:
        '''
        /home/workspace/jogrady/heQTL/software/bcftools-1.15.1/bcftools query -f '%ID\t%CHROM\t%POS\t%POS\t%INFO\n' {input.vcf} | sed 's/chr//g' > /home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/snp_density.txt

        Rscript {input.script} /home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/snp_density.txt {output.plot} {output.zoom1} {output.zoom2} 
        '''
