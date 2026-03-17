# def get_star_consensus_index_input(wc):
#     ref_name = config['datasets'][wc.cond]['reference_genotype']
#     parent2_accessions = config['datasets'][wc.cond]['parent2_accessions']
#     input_ = {
#         'fasta_fn': ancient(config['genome_fasta'][ref_name]),
#         'organellar_fasta_fn': ancient(config['organellar_fasta']),
#         'gtf_fn': ancient(config['genome_gff'][ref_name]),
#         'organellar_gtf_fn': ancient(config['organellar_gff']),
#     }
#     if len(parent2_accessions) > 1:
#         input_['vcf_fn'] = '../../annotations/diploid/vcf/{cond}.consensus.vcf'
#     return input_ 


# rule build_STAR_index_consensus:
#     '''Create the index required for alignment with STAR'''
#     input:
#         unpack(get_star_consensus_index_input)
#     output:
#         directory('../../annotatios/reference/{cond}_STAR_index')
#     threads: 12
#     resources:
#         mem_mb=12 * 1024,
#         queue='ioheavy'
#     params:
#         overhang=150,
#         vcf_flag=lambda wc, input: f'--genomeTransformVCF {input.vcf_fn}' if hasattr(input, 'vcf_fn') else '',
#         transform_flag=lambda wc, input: '--genomeTransformType Haploid' if hasattr(input, 'vcf_fn') else '',
#     conda:
#         'env_yml/star.yml'
#     shell:
#         '''
#         mkdir {output};
#         cat {input.fasta_fn} {input.organellar_fasta_fn} > \
#           ../../annotations/diploid/{wildcards.cond}_fullref.fasta
#         cat {input.gtf_fn} {input.organellar_gtf_fn} > \
#           ../../annotations/diploid/{wildcards.cond}_fullref.gff
#         STAR \
#           --runThreadN {threads} \
#           --runMode genomeGenerate \
#           --genomeDir {output} \
#           --genomeFastaFiles ../../annotations/diploid/{wildcards.cond}_fullref.fa \
#           --sjdbGTFfile ../../annotations/diploid/{wildcards.cond}_fullref.gff \
#           {params.vcf_flag} \
#           {params.transform_flag} \
#           --genomeSAindexNbases 12 \
#           --sjdbOverhang {params.overhang}

#         rm ../../annotations/diploid/{wildcards.cond}_fullref.fasta \
#            ../../annotations/diploid/{wildcards.cond}_fullref.gff
#         '''


# rule barcode_trim:
#     '''
#     Trim the cDNA sequence from the fwd read to yield only the barcode
#     since this is the format required by STARsolo 
    
#     This seems to be a bit of a waste of the cDNA info in read 1.
#     '''
#     input:
#         read='raw_data/{sample_name}.1.fastq.gz',
#     output:
#         read_barcode='trimmed_data/{sample_name}.1_barcode.fastq.gz',
#     threads: 4
#     resources:
#         mem_mb=2048,
#         queue='ioheavy'
#     conda:
#         'env_yamls/seqkit.yaml'
#     shell:
#         '''
#         seqkit subseq -r 1:28 {input.read} | gzip > {output.read_barcode}
#         '''


# def sample_name_subset(cond):
#     sample_names = glob_wildcards(
#         'raw_data/{sample_name}.1.fastq.gz'
#     ).sample_name
#     cond_sample_names = [sn for sn in sample_names if sn.rsplit('_', 1)[0] == cond]
#     return sorted(cond_sample_names)


# def expand_sample_name_from_cond(pattern):
#     def _expand_sn(wc):
#         return expand(
#             pattern,
#             sample_name=sample_name_subset(wc.cond)
#         )
#     return _expand_sn



# rule STARsolo_consensus:
#     '''
#     map reads with STAR spliced aligner
    
#     Read and mate are swapped in input specification compared to regular STAR run
#     as read 2 contains the cDNA and read 1 is the barcode.
#     '''
#     input:
#         read_barcode=expand_sample_name_from_cond('trimmed_data/{sample_name}.1_barcode.fastq.gz'),
#         mate=expand_sample_name_from_cond('raw_data/{sample_name}.2.fastq.gz'),
#         index='../../annotations/reference/{cond}_STAR_index',
#         barcode_whitelist=ancient(config['barcode_whitelist']),
#     output:
#         bam='aligned_data/reference/{cond}.sorted.bam',
#         bai='aligned_data/reference/{cond}.sorted.bam.bai',
#         stats='aligned_data/reference/{cond}.sorted.bamstats',
#         sjdb='aligned_data/reference/{cond}.sjdb.tsv',
#         solo_output=directory('aligned_data/reference/{cond}_starsolo')
#     params:
#         read_barcode=lambda wc, input: ','.join(f'${{TOPDIR}}/{fn}' for fn in input.read_barcode),
#         mate=lambda wc, input: ','.join(f'${{TOPDIR}}/{fn}' for fn in input.mate),
#         sort_mem=lambda wc, resources: (resources.mem_mb - 4096) * 1_000_000,
#         n_files=lambda wc, threads: threads * 150 + 200,
#         transform_flag=lambda wc: '--genomeTransformOutput SAM SJ Quant' if len(config['datasets'][wc.cond]['parent2_accessions']) > 1 else ''
#     log:
#         progress='logs/{cond}.STAR_progress.log',
#         final='logs/{cond}.STAR_final.log',
#         main='logs/{cond}.STAR.log'
#     threads: 24
#     resources:
#         mem_mb=lambda wildcards, threads: (threads + 4) * 2048,
#         queue='ioheavy'
#     conda:
#         'env_yamls/star.yaml'
#     shell:
#         '''
#         TOPDIR=$(pwd)
#         STAR_TMP_DIR="aligned_data/reference/{wildcards.cond}.tmpdir"
#         mkdir -p $STAR_TMP_DIR
#         cd $STAR_TMP_DIR
#         ulimit -n {params.n_files}
#         STAR \
#           --runThreadN {threads} \
#           --genomeDir "$TOPDIR/{input.index}" \
#           --readFilesIn "{params.mate}" "{params.read_barcode}" \
#           --readFilesCommand "zcat" \
#           --soloType "CB_UMI_Simple" \
#           --soloCBwhitelist "$TOPDIR/{input.barcode_whitelist}" \
#           --soloUMIlen 12 \
#           --soloUMIdedup "1MM_Directional_UMItools" \
#           --outFilterMultimapNmax 1 \
#           --outFilterIntronMotifs RemoveNoncanonical \
#           --alignSJoverhangMin 12 \
#           --alignSJDBoverhangMin 4 \
#           --outFilterMismatchNmax 4 \
#           --alignIntronMin 60 \
#           --alignIntronMax 20000 \
#           --outSAMunmapped Within \
#           --outSAMtype BAM SortedByCoordinate \
#           --outBAMsortingBinsN 150 \
#           --limitBAMsortRAM {params.sort_mem} \
#           --outSAMattributes NH HI AS nM CB UB UR sS sQ \
#           {params.transform_flag}

#         cd $TOPDIR
#         mv ${{STAR_TMP_DIR}}/Aligned.sortedByCoord.out.bam {output.bam}
#         mv ${{STAR_TMP_DIR}}/Log.progress.out {log.progress}
#         mv ${{STAR_TMP_DIR}}/Log.final.out {log.final}
#         mv ${{STAR_TMP_DIR}}/Log.out {log.main}
#         mv ${{STAR_TMP_DIR}}/SJ.out.tab {output.sjdb}
#         mv ${{STAR_TMP_DIR}}/Solo.out {output.solo_output}
          
#         samtools index {output.bam}
#         samtools flagstat {output.bam} > {output.stats}
#         rm -rf $STAR_TMP_DIR
#         '''


# rule cb_whitelist:
#     input:
#         solo_output='aligned_data/reference/{cond}_starsolo',
#     output:
#         whitelist='aligned_data/reference/{cond}.cb_whitelist.txt',
#     threads: 1
#     resources:
#         mem_mb=10_000
#     conda:
#         'env_yamls/pomegranate.yaml'
#     shell:
#         '''
#         python ../scripts/cb_whitelist.py \
#           -s {input.solo_output} \
#           -o {output.whitelist}
#         '''


rule cellsnp_lite:
    input:
        bam='../../aligned_data/reference/{cond}.sorted_Chr4.bam',
        barcodes='../../aligned_data/reference/{cond}.cb_whitelist.txt',
        vcf='../../annotations/diploid/vcf/{cond}.cellsnplite_fixed_Chr4.vcf'
    output:
        vcf='../../cb_snp_counts/{cond}/cellSNP.base.vcf',
        bc='../../cb_snp_counts/{cond}/cellSNP.samples.tsv',
        ad='../../cb_snp_counts/{cond}/cellSNP.tag.AD.mtx',
        dp='../../cb_snp_counts/{cond}/cellSNP.tag.DP.mtx',
    threads: 8
    resources:
        mem_mb=10_000
    conda:
        'env_yml/cellsnp.yml'
    shell:
        '''
        cellsnp-lite \
          -s {input.bam} \
          -b {input.barcodes} \
          -R {input.vcf} \
          -p {threads} \
          --minCount 0 \
          -O ../../cb_snp_counts/{wildcards.cond}
        '''
        
rule cellSNP_to_RTIGER:
    input:
        samples="../../cb_snp_counts/col0xdb1/cellSNP.samples.tsv",
        vcf="../../cb_snp_counts/col0xdb1/cellSNP.base.vcf",
        ad="../../cb_snp_counts/col0xdb1/cellSNP.AD.mtx",
        dp="../../cb_snp_counts/col0xdb1/cellSNP.DP.mtx"
    output:
        directory("../../RTIGER/{cond}/samples")
    threads: 1
    resources:
    	mem_mb=4096
    conda:
    	'env_yml/rtiger.yml'
    script:
    	"../scripts/rtiger_input.R"

rule run_rtiger:
	input:
		"../../RTIGER/{cond}/samples/AAACCCACACACCTTC.tsv"
	output:
		"../../RTIGER/{cond}/myres_object.rds"
	threads: 4
	resources:
		mem_mb=8_000
	conda:
		"env_yml/rtiger.yml"
	script:
		"../scripts/rtiger.R"

rule assign_genotypes:
    input:
        vcf='../../annotations/diploid/vcf/{cond}.cellsnplite.vcf.gz',
        bc='../../cb_snp_counts/{cond}/cellSNP.samples.tsv',
        ad='../../cb_snp_counts/{cond}/cellSNP.tag.AD.mtx',
        dp='../../cb_snp_counts/{cond}/cellSNP.tag.DP.mtx',
    output:
        geno='../../cb_genotypes/{cond}.geno.tsv'
    threads: 8
    resources:
        mem_mb=30_000
    conda:
        'env_yml/pysam.yml'
    shell:
        '''
        python ../scripts/assign_genotypes.py -p {threads} \
          -c cb_snp_counts/{wildcards.cond}/ \
          -v {input.vcf} \
          -o {output.geno}
        '''


# rule demux_genotypes:
#     input:
#         bam='aligned_data/reference/{cond}.sorted.bam',
#         geno='cb_genotypes/{cond}.geno.tsv',
#     output:
#         read1='demuxed_data/{cond}.{parent1}__vs__{parent2}.1_barcode.fastq.gz',
#         read2='demuxed_data/{cond}.{parent1}__vs__{parent2}.2.fastq.gz',
#     threads: 9
#     resources:
#         mem_mb=10_000,
#         queue='ioheavy'
#     params:
#         threads=lambda wc, threads: threads // 3,
#     conda:
#         'env_yamls/samtools.yaml'
#     shell:
#         '''
#         awk '$2 == "{wildcards.parent2}" {{print $1}}' {input.geno} > "demuxed_data/{wildcards.parent2}.cb.tsv"
#         samtools view -@ {params.threads} -hb -D "CB:demuxed_data/{wildcards.parent2}.cb.tsv" {input.bam} | \
#         samtools sort -n -m 1G -@ {params.threads} -o - - | \
#         samtools fastq -@ {params.threads} \
#           -0 {output.read2} \
#           --i1 {output.read1} \
#           --barcode-tag sS \
#           --quality-tag sQ \
#           --index-format "i*" 
#         rm "demuxed_data/{wildcards.parent2}.cb.tsv"
#         '''
