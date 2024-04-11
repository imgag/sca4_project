wildcard_constraints:
    sample = "[A-Za-z0-9]+" 

rule consensusreadset:
    input:
        lambda w: samples_pb[w.sample]
    output:
        "pacbio/{sample}.ccs.xml"
    conda:
        "env/pbtools.yml"
    shell:
        """
        dataset merge {output} {input} --name {wildcards.sample}
        """

rule alignment:
    input:
        ref = lambda w: ref[w.ref],
        xml = "pacbio/{sample}.ccs.xml"
    output:
        bam = "pacbio/alignment_{ref}/{sample}/{sample}.{ref}.bam"
    conda:
        "env/pbmm2.yml"
    log:
        "logs/{sample}_{ref}_pbmm2.log"
    threads:
        16
    shell:
        """
        pbmm2 align --preset CCS --unmapped -j {threads} {input.ref} {input.xml} {output} --sort -j4 -J 2 >{log} 2>&1
        """

rule qualimap:
    input:
        "pacbio/alignment_{ref}/{sample}/{sample}.{ref}.bam"
    output:
        'pacbio/alignment_{ref}/{sample}/qualimap_{sample}/genome_results.txt'
    log:
        "logs/{sample}_{ref}_qualimap.log"
    threads:
        8
    conda:
        "env/qualimap.yml"
    shell:
        """
        qualimap bamqc \
        -bam {input} \
        --paint-chromosome-limits \
        -nt {threads} \
        -outdir $(dirname {output}) \
        --java-mem-size=24G \
        """
    