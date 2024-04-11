rule bam_to_fastq:
    input:
        "pacbio/{sample}.ccs.xml"
    output:
        "pacbio/{sample}.ccs.fastq.gz"
    conda:
        "env/pbtools.yml"
    shell:
        """
        bam2fastq -o pacbio/{wildcards.sample}.ccs {input}
        """

rule jellyfish:
    input:
        "pacbio/{sample}.ccs.fastq.gz"
    output:
        "pacbio/genomescope/{sample}/reads.jf"
    conda:
        "env/genomescope.yml"
    threads:
        6
    shell:
        "jellyfish count -C -m 21 -s 1000000000 -t {threads} <(zcat {input}) -o {output}"

rule jellyfish_histo:
    input:
        "pacbio/genomescope/{sample}/reads.jf"
    output:
        "pacbio/genomescope/{sample}/reads.histo"
    conda:
        "env/genomescope.yml"
    threads: 
        6
    shell:
        "jellyfish histo -t{threads} {input} > {output}"

rule genomescope:
    input:
        "pacbio/genomescope/{sample}/reads.histo"
    output:
        "pacbio/genomescope/{sample}/summary.txt"
    conda:
        "env/genomescope.yml"
    shell:
        """
        genomescope2 \
            -i {input} \
            -o $(dirname {output}) \
            -p 2 \
            -k 21
        """

rule hifiasm:
    input:
        "pacbio/{sample}.ccs.fastq.gz"
    output:
        "pacbio/assembly_hifiasm/{sample}/{sample}.asm_hifi.bp.p_ctg.gfa"
    conda:
        "env/hifiasm.yml"
    log: 
        "logs/{sample}_hifiasm.log"
    threads:
        28
    shell:
        """
        hifiasm  \
            -o pacbio/assembly_hifiasm/{wildcards.sample}/{wildcards.sample}.asm_hifi \
            -t {threads} \
            {input} > {log} 2>&1
        """

rule create_fasta:
    input:
        "pacbio/assembly_hifiasm/{sample}/{sample}.{type}.gfa"
    output:
        "pacbio/assembly_hifiasm/{sample}/{sample}.{type}.fasta"
    threads:
        1
    shell:
        """
        awk '/^S/{{print ">"$2;print $3}}' {input} > {output}
        """

rule quast:
    input:
        fa = "pacbio/assembly_hifiasm/{sample}/{sample}.asm_hifi.bp.p_ctg.fasta",
        fq = "pacbio/{sample}.ccs.fastq.gz"
    output:
        directory("pacbio/assembly_hifiasm/{sample}/quast")
    conda:
        "env/quast.yml"
    log:
        "logs/{sample}_quast.log"
    threads:
        8
    shell:
        """
        quast.py -b \
            --circos \
            --eukaryote \
            --glimmer \
            --threads {threads} \
            -o {output} \
            {input.fa} > {log} 2>&1
        """
        