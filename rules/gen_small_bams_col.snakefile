rule chr4_fasta:
    input:
        db='../../annotations/Db-1.fasta.gz',
        cc='../../annotations/Col-CC.fasta.gz'
    output:
        db='../../annotations/Db-1_Chr4.fasta',
        cc='../../annotations/Col-CC_Chr4.fasta'
    conda:
        'env_yml/pysam.yml'
    shell:
        '''
        zcat {input.db} | awk '/^>/ {{flag = ($0 ~ /^>Chr4/)}} flag' > {output.db}
        zcat {input.cc} | awk '/^>/ {{flag = ($0 ~ /^>Chr4/)}} flag' > {output.cc}
        '''

rule chr4_bam:
    input:
        '../../annotations/col0xdb1.sorted.bam'
    output:
        bam='../../annotations/col0xdb1_Chr4.sorted.bam',
        bai='../../annotations/col0xdb1_Chr4.sorted.bam.bai'
    threads:
        4
    conda:
        'env_yml/pysam.yml'
    shell:
        '''
        samtools view -b -@ {threads} {input} chr4 > {output.bam}
        samtools index {output.bam}
        '''

rule subsample_bam:
    input:
        '../../annotations/col0xdb1.sorted.bam'
    output:
        bam='../../annotations/col0xdb1_subsample.sorted.bam',
        bai='../../annotations/col0xdb1_subsample.sorted.bam.bai'
    threads:
        4
    conda:
        'env_yml/pysam.yml'
    shell:
        '''
        samtools view -bs 123.2 -@ {threads} {input} > {output.bam}
        samtools index {output.bam}
        '''

rule chr4_gff:
    input:
        db='../../annotations/Db-1.gff',
        cc='../../annotations/Col-CC.gff'
    output:
        db='../../annotations/Db-1_Chr4.gff',
        cc='../../annotations/Col-CC_Chr4.gff'
    conda:
        'env_yml/pysam.yml'
    shell:
        '''
        cat {input.db} | grep "^Chr4" > {output.db}
        cat {input.cc} | grep "^Chr4" > {output.cc}
        '''
