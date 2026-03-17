def get_wga_input(wc):
    parent1, parent2 = wc.comp.split('__vs__')
    return {
        'parent1': ancient(config['genome_fasta'][parent1]),
        'parent2': ancient(config['genome_fasta'][parent2]),
    }


rule minimap2_wga:
    input:
        unpack(get_wga_input)
    output:
        '../../annotations/diploid/wga/{comp}.bam'
    threads:
        8
    resources:
        mem_mb=lambda wc, threads: 2048 * threads
    conda:
        'env_yml/minimap.yml'
    shell:
        '''
        minimap2 --eqx -t {threads} -ax asm20 -z1000,100 --no-end-flt \
           {input.parent1} {input.parent2}  | \
        samtools view -bS | \
        samtools sort -o - - > {output}
        samtools index {output}
        '''


def get_syri_input(wc):
    input_ = get_wga_input(wc)
    input_['bam'] = '../../annotations/diploid/wga/{comp}.bam'
    return input_


rule run_syri:
    input:
        unpack(get_syri_input)
    output:
        vcf='../../annotations/diploid/vcf/{comp}.syri.vcf',
        syri='../../annotations/diploid/vcf/{comp}.syri.out',
    threads:
        1
    resources:
        mem_mb=8192
    conda:
        'env_yml/msyd.yml'
    params:
        acc=lambda wc: wc.comp.split('__vs__')[1]
    shell:
        '''
        syri -F B -f --hdrseq --dir ../../annotations/diploid/vcf --prefix {wildcards.comp}. \
          -c {input.bam} -q {input.parent2} -r {input.parent1} --samplename {params.acc}
        '''

rule plotsr_syri:
    input:
        syri='../../annotations/diploid/vcf/{comp}.syri.out',
        genomes='../../annotations/genomes.txt',
    output:
        '../../results/{comp}_syri_plot.png'
    threads:
        1
    resources:
        mem_mb=8192
    conda:
        'env_yml/msyd.yml'
    params:
        acc=lambda wc: wc.comp.split('__vs__')[1]
    shell:
        '''
        plotsr --sr {input.syri} --genomes {input.genomes} -o {output}
        '''

def get_msyd_input(wc):
    accs = config['datasets'][wc.cond]['parent2_accessions']
    ref_name = config['datasets'][wc.cond]['reference_genotype']
    return {
        'bams': expand('../../annotations/diploid/wga/{ref_name}__vs__{accession}.bam', ref_name=ref_name, accession=accs),
        'syri': expand('../../annotations/diploid/vcf/{ref_name}__vs__{accession}.syri.out', ref_name=ref_name, accession=accs),
        'vcf': expand('../../annotations/diploid/vcf/{ref_name}__vs__{accession}.syri.vcf', ref_name=ref_name, accession=accs),
        'fasta': [config['genome_fasta'][acc] for acc in accs],
    }


rule msyd_input:
    input:
        unpack(get_msyd_input)
    output:
        cfg=temp('../../annotations/diploid/vcf/{cond}.msyd_config.tsv')
    run:
        ref_name = config['datasets'][wildcards.cond]['reference_genotype']
        with open(output.cfg, 'w') as f:
            f.write('#name\taln\tsyri\tvcf\tgenome\n')
            for accession in config['datasets'][wildcards.cond]['parent2_accessions']:
                f.write(
                    f'{accession}\t'
                    f'../../annotations/diploid/wga/{ref_name}__vs__{accession}.bam\t'
                    f'../../annotations/diploid/vcf/{ref_name}__vs__{accession}.syri.out\t'
                    f'../../annotations/diploid/vcf/{ref_name}__vs__{accession}.syri.vcf\t'
                    f'{config["genome_fasta"][accession]}\n'
                )


rule run_msyd:
    input:
        cfg='../../annotations/diploid/vcf/{cond}.msyd_config.tsv',
        ref=lambda wc: ancient(config['genome_fasta'][config['datasets'][wc.cond]['reference_genotype']])
    output:
        pff=r'../../annotations/diploid/vcf/{cond,\w+}.pansyn.pff',
        vcf=r'../../annotations/diploid/vcf/{cond,\w+}.vcf',
    conda:
        'env_yml/msyd.yml'
    shell:
        '''
        msyd call --core -i {input.cfg} -r {input.ref} -o {output.pff} -m {output.vcf}
        '''


rule convert_vcf_for_star:
    input:
        '../../annotations/diploid/vcf/{comp}.syri.vcf'
    output:
        '../../annotations/diploid/vcf/{comp}.snvs.vcf'
    conda:
        'env_yml/pysam.yml'
    resources:
        mem_mb=8_000
    shell:
        '''
        python ../scripts/convert_syri_output.py {input} {output}
        '''


rule filter_snps_for_star_consensus:
    input:
        vcf=r'../../annotations/diploid/vcf/{cond}.vcf'
    output:
        vcf=r'../../annotations/diploid/vcf/{cond,\w+}.consensus.vcf',
    conda:
        'env_yml/pysam.yml'
    params:
        min_ac=lambda wc: len(config['datasets'][wc.cond]['parent2_accessions']) // 2 + 1,
        max_indel_size=50,
    shell:
        r'''
        bcftools annotate \
          --exclude 'ALT ~ "CORESYN"' \
          --remove "FORMAT/CHR,FORMAT/START,FORMAT/END,INFO/PID" \
          {input.vcf} |
        grep -v "^##ALT" | grep -v "^##INFO" | \
        bcftools sort | \
        bcftools filter -S0 -e 'GT=="."' | \
        bcftools view -G -M2 \
          --min-ac {params.min_ac} \
          -e "STRLEN(REF)>{params.max_indel_size} || STRLEN(ALT)>{params.max_indel_size}" \
        > {output.vcf}
        '''


rule filter_snps_for_cellsnplite:
    input:
        vcf='../../annotations/diploid/vcf/{cond}.vcf'
    output:
        vcf=r'../../annotations/diploid/vcf/{cond,\w+}.cellsnplite.vcf.gz'
    conda:
        'env_yml/pysam.yml'
    shell:
        r'''
        bcftools annotate \
          --exclude 'ALT ~ "CORESYN"' \
          --remove "FORMAT/CHR,FORMAT/START,FORMAT/END,INFO/PID" \
          {input.vcf} |
        grep -v "^##ALT" | grep -v "^##INFO" | \
        bcftools sort | \
        bcftools filter -S0 -e 'GT=="."' | \
        bcftools view -M2 \
          -i 'TYPE="snp"' \
          -O z \
        > {output.vcf}
        tabix -p vcf {output.vcf}
        '''
