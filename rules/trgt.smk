# binary executable locations
pbmm2   = '/mnt/storage2/users/ahgrosc1/environments/envs/pb_tools/bin/pbmm2'
trgt    = '/mnt/storage2/users/ahgrosc1/tools/bin/trgt-v0.4.0-linux_x86_64'
trvz    = '/mnt/storage2/users/ahgrosc1/tools/bin/trvz-v0.4.0-linux_x86_64'


rule trgt:
    input:
        genome = ref['GRCh38'],
        repeats = lambda wc: resources[wc.set],
        reads = "pacbio/alignment_GRCh38/{sample}/{sample}.GRCh38.bam",
    output:
        vcf = "pacbio/repeats_trgt/{sample}.{set}.GRCh38.vcf.gz",
        bam = "pacbio/repeats_trgt/{sample}.{set}.GRCh38.spanning.bam"
    params:
        trgt = trgt
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        {params.trgt} \
            --genome {input.genome} \
            --repeats {input.repeats} \
            --reads {input.reads} \
            --output-prefix pacbio/repeats_trgt/{wildcards.sample}.{wildcards.set}.GRCh38
        samtools index {output.bam}
        """
rule index_spanning:
    input:
        bam = "pacbio/repeats_trgt/{sample}.{set}.GRCh38.spanning.bam"
    output:
        bai = "pacbio/repeats_trgt/{sample}.{set}.GRCh38.spanning.bam.bai"
    shell:
        """
        samtools index {input.bam}
        """

checkpoint find_expanded:
    input:
        vcf = "pacbio/repeats_trgt/{sample}.{set}.GRCh38.vcf.gz"
    output:
        tsv = "pacbio/repeats_trgt/{sample}.{set}.expanded_loci.tsv",
        folder = directory("pacbio/repeats_trgt/{sample}.{set}.plot") 
    conda:
        "env/env_r.yml"
    script:
        "find_expanded.R"

rule trvz:
    input:
        genome = ref['GRCh38'],
        repeats = lambda wc: resources[wc.set],
        vcf = "pacbio/repeats_trgt/{sample}.{set}.GRCh38.vcf.gz",
        bam = "pacbio/repeats_trgt/{sample}.{set}.GRCh38.spanning.bam",
        bai = "pacbio/repeats_trgt/{sample}.{set}.GRCh38.spanning.bam.bai"
    output:
        image = "pacbio/repeats_trgt/{sample}.{set}.plot/{id}/{id}.{type}.{vis}.{ext}"
    params:
        trvz = trvz
    shell:
        """
        {params.trvz} \
            --genome {input.genome} \
            --repeats {input.repeats} \
            --vcf {input.vcf} \
            --plot-type {wildcards.type} \
            --spanning-reads {input.bam} \
            --repeat-id {wildcards.id} \
            --show {wildcards.vis} \
            --image {output.image} 
        """

def input_repeat_ids(wildcards):
    checkpoint_output = checkpoints.find_expanded.get(**wildcards).output.folder
    #print(checkpoint_output)
    repeat_folders = glob_wildcards(os.path.join(checkpoint_output, "{f}", "repeat_locus_info.tsv")).f
    #print(repeat_folders)
    return expand("pacbio/repeats_trgt/{{sample}}.{{set}}.plot/{id}/{id}.{{type}}.{{vis}}.{{ext}}",
        id = repeat_folders)
    
rule plot_aggregate:
    input: 
        input_repeat_ids
    output:
        "pacbio/repeats_trgt/plot.{sample}.{set}.{type}.{vis}.{ext}.done"
    shell:
        """
        touch {output}
        """