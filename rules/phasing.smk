rule longphase_phase:
    input:
        vcf_snp = "pacbio/variant_calling/{sample}_{ref}_deepvariant/{sample}.dv.vcf.gz",
        bam = "pacbio/alignment_{ref}/{sample}/{sample}.{ref}.bam",
        ref = lambda wc: ref[wc.ref]
    output:
        vcf = "pacbio/phasing/{sample}_{ref}/{sample}.phased.vcf.gz",
        vcf_sv = "pacbio/phasing/{sample}_{ref}/{sample}.phased_SV.vcf.gz"
    threads:
        4
    log:
        "logs/{sample}_longphase_phase_{ref}.log"
    params:
        longphase = "rules/env/longphase_linux-x64",
        sv = lambda wc: "--sv-file pacbio/variant_calling/{sample}_{ref}_dipdiff/hapdiff_unphased.vcf.gz".format(sample = wc.sample, ref = wc.ref) if os.path.isfile("pacbio/variant_calling/{sample}_{ref}_dipdiff/hapdiff_unphased.vcf.gz".format(sample = wc.sample, ref = wc.ref)) else ""
    shell:
        """
        {params.longphase} phase \
        --snp-file {input.vcf_snp} \
        {params.sv} \
        --bam-file {input.bam} \
        --pb \
        --reference {input.ref} \
        --threads {threads} \
        --out-prefix pacbio/phasing/{wildcards.sample}_{wildcards.ref}/{wildcards.sample}.phased \
        > {log} 2>&1

        bgzip pacbio/phasing/{wildcards.sample}_{wildcards.ref}/{wildcards.sample}.phased.vcf
        bgzip pacbio/phasing/{wildcards.sample}_{wildcards.ref}/{wildcards.sample}.phased_SV.vcf
        tabix {output.vcf}
        tabix {output.vcf_sv}
        """

rule longphase_haplotag:
    input:
        vcf_snp = "pacbio/phasing/{sample}_{ref}/{sample}.phased.vcf.gz",
        bam = "pacbio/alignment_{ref}/{sample}/{sample}.{ref}.bam",
    output:
        "pacbio/alignment_{ref}/{sample}/{sample}.{ref}.haplotagged.bam"
    threads:
        6
    log:
        "logs/{sample}_{ref}_longphase_haplotagging.log"
    params:
        longphase = "rules/env/longphase_linux-x64"
    shell:
        """
        {params.longphase} haplotag \
        --snp-file {input.vcf_snp} \
        --bam-file {input.bam} \
        --output-threads {threads} \
        --out-prefix pacbio/alignment_{wildcards.ref}/{wildcards.sample}/{wildcards.sample}.{wildcards.ref}.haplotagged \
        > {log} 2>&1

        samtools index {output}
        """

rule fix_vcf_name:
    input:
        "pacbio/phasing/{sample}_{ref}/{sample}.phased.vcf.gz",
    output:
        "pacbio/phasing/{sample}_{ref}/{sample}.phased.renamed.vcf.gz"
    conda:
        "env/glnexus.yml"
    shell:
        """
        bcftools reheader \
            --sample <(echo {wildcards.sample}) \
            {input} \
        | bcftools view \
            --output-type z \
            --output {output} \
            -
        tabix {output}
        """

rule phasing_stats:
    input:
        vcf = "pacbio/phasing/{sample}_{ref}/{sample}.phased.renamed.vcf.gz",
    output:
        tsv = "pacbio/phasing/{sample}_{ref}/{sample}_{ref}.phasestats.tsv",
        blocks = "pacbio/phasing/{sample}_{ref}/{sample}_{ref}.phaseblocks.tsv",
        gtf = "pacbio/phasing/{sample}_{ref}/{sample}_{ref}.phaseblocks.gtf"
    conda:
        "env/whatshap.yml"
    log:
        "logs/{sample}_{ref}_whatshap_stats.log"
    shell:
        """
        whatshap stats {input.vcf} \
            --tsv={output.tsv} \
            --block-list={output.blocks} \
            --gtf={output.gtf} > {log} 2>&1
        """