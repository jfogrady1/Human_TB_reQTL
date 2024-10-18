#######################################################################################
####### snakemake script for checking for amplification bias in sample 3 at timepoint 1 (SRR12609781)
##################################################################################################
#rule all:
 #    input:
  #    "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/DV_snp_tab.h5",
   #   "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/intersecting_snps/SRR12609781.keep.bam",
    #  '/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/map2/SRR12609781.sort.bam',
     # "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/filter_remapped/SRR12609781.keep.bam",
      #"/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/merge/SRR12609781.keep.merge.sort.bam",
      #"/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/dup_filtered/sorted/SRR12609781.keep.merge.sort.no_dups.bam",
      #"/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/dup_filtered/sorted/SRR12609781.mbv_output.txt"

autosome = [str(i) for i in range(1, 23)]

ruleorder: vcf2hd5 > find_intersecting_snps > Alignment_second_round > filter_remapped_reads
#rule 1

rule vcf2hd5:
     input: 
        chrominfo ="/home/workspace/jogrady/heQTL/data/ref_genome/GRCh38.chrominfo.txt", 
        vcf = "/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf.gz",
        script = "/home/workspace/jogrady/heQTL/software/WASP-master/snp2h5/snp2h5"
     output: 
        hd5 = "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/DV_snp_tab.h5",
        snp_index = "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/DV_snp_index.h5",
        haplotype = "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/DV_haplotypes.h5"
        
     shell:
        '''
         {input.script} --chrom {input.chrominfo} \
         --format vcf \
         --haplotype /home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/DV_haplotypes.h5 \
         --snp_index /home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/DV_snp_index.h5 \
         --snp_tab {output.hd5} \
         {input.vcf}
        '''
        
# we have doen rule 2 (mapping of reads)



# Rule3
rule find_intersecting_snps:
  input:
    bam = '/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/SRR12609781.bam',
    haplotype = '/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/DV_haplotypes.h5',
    index = '/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/DV_snp_index.h5',
    snp_tab = '/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/DV_snp_tab.h5',
    samples = '/home/workspace/jogrady/heQTL/data/ref_genome/WASP_my_samples.txt'

  output:
    bam_keep = "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/intersecting_snps/SRR12609781.keep.bam",
    reads = expand(f'/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/intersecting_snps/SRR12609781.remap.fq{{N}}.gz', N=(1,2))
  
  shell:
    """
    python /home/workspace/jogrady/heQTL/software/WASP-master/mapping/find_intersecting_snps.py \
              --is_paired_end \
              --is_sorted \
              --output_dir find_intersecting_snps \
              --snp_tab {input.snp_tab} \
              --snp_index {input.index} \
              --haplotype {input.haplotype} \
              --samples {input.samples} \
              --output_dir /home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/intersecting_snps/ \
              {input.bam}
    """

# Rule 4 - realignment
rule Alignment_second_round:
    input:
        genome = '/home/workspace/jogrady/heQTL/work/RNA_seq/STAR_DIR_updated',
        reads = expand(f'/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/intersecting_snps/SRR12609781.remap.fq{{N}}.gz', N=(1,2)),

    params:
        prefix = 'SRR12609781'
    output:
        aligned = '/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/map2/SRR12609781_Aligned.sortedByCoord.out.bam',
        finallog = '/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/map2/SRR12609781_Log.final.out',
        interlog = '/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/map2/SRR12609781_Log.progress.out',
        initiallog = '/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/map2/SRR12609781_Log.out',
        sorted = '/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/map2/SRR12609781.sort.bam'

    threads: 20

    shell:
        """
        /home/workspace/jogrady/heQTL/software/STAR-2.7.10b/bin/Linux_x86_64/STAR --genomeLoad NoSharedMemory \
        --genomeDir {input.genome} --runThreadN {threads} \
        --readFilesIn {input.reads[0]} {input.reads[1]} \
        --readFilesCommand gunzip -c \
        --outFileNamePrefix /home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/map2/{params.prefix}_ \
        --outSAMtype BAM SortedByCoordinate \
        --limitBAMsortRAM 5000000000 
        
        /home/workspace/jogrady/heQTL/software/samtools-1.15.1/samtools view -b -q 10 /home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/map2/SRR12609781_Aligned.sortedByCoord.out.bam > /home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/map2/SRR12609781.bam

        /home/workspace/jogrady/heQTL/software/samtools-1.15.1/samtools sort -o /home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/map2/SRR12609781.sort.bam /home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/map2/SRR12609781.bam

        /home/workspace/jogrady/heQTL/software/samtools-1.15.1/samtools index /home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/map2/SRR12609781.sort.bam
        """
        
rule filter_remapped_reads:
  input:
    remapped_bam = '/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/map2/SRR12609781.sort.bam'

  output:
    bam_keep = "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/filter_remapped/SRR12609781.keep.bam"
  
  shell:
    """
    /home/workspace/jogrady/heQTL/software/samtools-1.15.1/samtools index {input.remapped_bam}
    
    python /home/workspace/jogrady/heQTL/software/WASP-master/mapping/filter_remapped_reads.py \
       /home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/intersecting_snps/SRR12609781.to.remap.bam \
       /home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/map2/SRR12609781.sort.bam \
       {output.bam_keep}
    """

rule form_mapp_filtered_aligned_read:
  input:
    remapped_bam_keep = "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/filter_remapped/SRR12609781.keep.bam",
    intersecting_bam_keep = "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/intersecting_snps/SRR12609781.keep.bam",
    samtools = "/home/workspace/jogrady/heQTL/software/samtools-1.15.1/samtools"
  output:
    merged_reads = "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/merge/SRR12609781.keep.merge.bam",
    sorted_merged_reads = "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/merge/SRR12609781.keep.merge.sort.bam"

  shell:
    """
    {input.samtools} merge {output.merged_reads} \
        {input.remapped_bam_keep}  \
        {input.intersecting_bam_keep}
    {input.samtools} sort -o  {output.sorted_merged_reads} \
              {output.merged_reads} 
    {input.samtools} index {output.sorted_merged_reads}
    """

rule filter_duplicates:
  input:
    script = "/home/workspace/jogrady/heQTL/software/WASP-master/mapping/rmdup_pe.py",
    sorted_bam = "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/merge/SRR12609781.keep.merge.sort.bam"

  output:
    filtered_bam = "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/dup_filtered/SRR12609781.keep.merge.sort.no_dups.bam"

  shell:
    """
    python {input.script} {input.sorted_bam} {output.filtered_bam}
    """

rule index_bam: 
  input:
    filtered_bam = "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/dup_filtered/SRR12609781.keep.merge.sort.no_dups.bam",
    samtools = "/home/workspace/jogrady/heQTL/software/samtools-1.15.1/samtools"

  output:
    sorted = "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/dup_filtered/sorted/SRR12609781.keep.merge.sort.no_dups.bam",
    sorted_index = "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/dup_filtered/sorted/SRR12609781.keep.merge.sort.no_dups.bam.bai"

  shell:
    """
    {input.samtools} sort -o  {output.sorted} \
              {input.filtered_bam} 
    {input.samtools} index {output.sorted}
    """


rule sample_check_1_WASP:
    input:
        vcf = multiext("/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered", ".vcf.gz", ".vcf.gz.tbi"),
        bam ="/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/dup_filtered/sorted/SRR12609781.keep.merge.sort.no_dups.bam",
        bai = "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/dup_filtered/sorted/SRR12609781.keep.merge.sort.no_dups.bam.bai"
    output:
        check = "/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/dup_filtered/sorted/SRR12609781.mbv_output.txt"
    singularity: "docker://jogrady/qtltools:1.3.1"

    shell:
        '''
        QTLtools mbv --vcf {input.vcf[0]} --bam {input.bam} --filter-mapping-quality 150 --filter-base-quality 20  --out {output.check} 
        '''