#!/usr/bin/env Rscript

# sshfs -o ro \
# bcc3:/Volumes/archive/deardenlab/HTS\ raw\ sequencing\ reads/phase_genomics \
# data/scaffolded

hyp = 'data/scaffolded/m_hyperodae/PGA_assembly.fasta'
aeth = 'data/scaffolded/m_aethiopoides/PGA_assembly.fasta'

biopython = 'shub://TomHarrop/py-containers:biopython_1.73'
minimap = 'shub://TomHarrop/align-utils:minimap2_2.17r941'
r = 'shub://TomHarrop/r-containers:r_4.0.0'
ragtag = 'shub://TomHarrop/assembly-utils:ragtag_1.0.1'
samtools = 'shub://TomHarrop/align-utils:samtools_1.10'

rule target:
    input:
        'wga/030_wga/wga.pdf'


rule plot_wga:
    input:
        query_fai = 'wga/020_oriented/m_hyperodae.fa.fai',
        ref_fai = 'wga/tmp/m_aethiopoides.fa.fai',
        paf = 'wga/030_wga/wga.paf'
    output:
        plot = 'wga/030_wga/wga.pdf'
    log:
        'wga/logs/plot_wga.log'
    singularity:
        r
    script:
        'src/plot_wga.R'


rule wga:
    input:
        fa = 'wga/020_oriented/m_hyperodae.fa',
        ref = 'wga/tmp/m_aethiopoides.mmi'
    output:
        pipe('wga/030_wga/wga.sam')
    log:
        'wga/logs/wga.log'
    threads:
        min(workflow.cores, 64) - 1
    singularity:
        minimap
    shell:
        'minimap2 '
        '-t {threads} '
        '-ax asm5 '
        '--MD '
        '{input.ref} '
        '{input.fa} '
        '> {output} '
        '2> {log}'


rule orient_scaffolds:
    input:
        fa = 'wga/tmp/m_hyperodae.fa',
        agp = 'wga/010_ragtag/ragtag.scaffolds.agp'
    output:
        fa = 'wga/020_oriented/m_hyperodae.fa'
    log:
        'wga/logs/orient_scaffolds.log'
    singularity:
        biopython
    script:
        'src/orient_scaffolds.py'

rule ragtag:
    input:
        'wga/tmp/m_hyperodae.fa.fai',
        'wga/tmp/m_aethiopoides.fa.fai',
        ref = 'wga/tmp/m_aethiopoides.fa',
        query = 'wga/tmp/m_hyperodae.fa'
    output:
        'wga/010_ragtag/ragtag.scaffolds.fasta',
        'wga/010_ragtag/ragtag.scaffolds.agp'
    params:
        wd = 'wga/010_ragtag'
    log:
        'wga/logs/ragtag.log'
    threads:
        min(workflow.cores, 64)
    container:
        ragtag
    shell:
        'ragtag.py scaffold '
        '-o {params.wd} '
        '-w '
        # '-r -g 101 '    # only add gaps 101 Ns or longer DOESN'T WORK
        '-t {threads} '
        '{input.ref} '
        '{input.query} '
        '&> {log}'



rule minimap_ref:
    input:
        'wga/tmp/m_aethiopoides.fa'
    output:
        'wga/tmp/m_aethiopoides.mmi'
    log:
        'wga/logs/prepare_ref.log'
    threads:
        3
    singularity:
        minimap
    shell:
        'minimap2 '
        '-x map-ont '
        '-d {output} '
        '{input} '
        '2> {log}'


rule tmp_ref:
    input:
        'data/scaffolded/{spec}/PGA_assembly.fasta'
    output:
        'wga/tmp/{spec}.fa'
    shell:
        'cp {input} {output}'


rule index_fa:
    input:
        '{path}/{file}.{ext}'
    output:
        '{path}/{file}.{ext}.fai'
    wildcard_constraints:
        ext = 'fasta|fa|fna'
    singularity:
        samtools
    shell:
        'samtools faidx {input}'


rule sam_to_paf:
    input:
        '{folder}/{indiv}.sam'
    output:
        '{folder}/{indiv}.paf'
    log:
        'wga/logs/sam_to_paf.{folder}.{indiv}.log'
    singularity:
        minimap
    shell:
        'paftools.js sam2paf '
        '{input} '
        '>{output} '
        '2>{log}'

rule sam_to_bam:
    input:
        '{folder}/{indiv}.sam'
    output:
        pipe('{folder}/{indiv}.bam')
    log:
        'wga/logs/sam_to_bam.{folder}.{indiv}.log'
    threads:
        1
    singularity:
        samtools
    shell:
        'samtools view -bh -u {input} '
        '>> {output} '
        '2> {log}'
