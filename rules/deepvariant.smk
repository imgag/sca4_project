rule snp_deepvariant:
  input: 
    bam = "pacbio/alignment_{ref}/{sample}/{sample}.{ref}.bam",
    ref = lambda w: ref[w.ref]
  output: 
    vcf = "pacbio/variant_calling/{sample}_{ref}_deepvariant/{sample}.dv.vcf.gz",
    gvcf = "pacbio/variant_calling/{sample}_{ref}_deepvariant/{sample}.dv.gvcf.gz"
  threads: 
    32
  log: 
    "logs/{sample}_{ref}_dv.log"
  resources:
    queue="research_srv023"
  shell:
    """
    set +u;
    docker run  \
      -v "$(realpath $(dirname {input.ref}))":"/mnt/ref" \
      -v "$(realpath $(dirname {input.bam}))":"/mnt/bam" \
      -v "$(realpath $(dirname {output.vcf}))":"/mnt/output" \
      -u `id -u`:`id -g` \
      google/deepvariant:1.4.0 \
      /opt/deepvariant/bin/run_deepvariant \
      --model_type=PACBIO \
      --ref=/mnt/ref/$(basename {input.ref}) \
      --reads=/mnt/bam/$(basename {input.bam}) \
      --output_vcf=/mnt/output/{wildcards.sample}.dv.vcf.gz \
      --output_gvcf=/mnt/output/{wildcards.sample}.dv.gvcf.gz \
      --num_shards={threads} \
      > {log} 2>&1
    """

rule annotate_megsap:
  input:  
    vcf = "pacbio/variant_calling/{sample}_{ref}_deepvariant/{sample}.dv.vcf.gz"
  output: 
    vcf = "pacbio/variant_calling/{sample}_{ref}_deepvariant/{sample}_var_annotated.vcf.gz"
  params:
    megsap_annotate = "php /mnt/storage2/megSAP/pipeline/src/Pipelines/annotate.php",
    ini = "annotations/hifi_grch38.ini"
  log:
    "logs/{sample}_{ref}_annotate.log"
  threads:
    4
  shell:
    """
    {params.megsap_annotate} \
      -out_name {wildcards.sample} \
      -out_folder $(dirname {output}) \
      -vcf {input.vcf} \
      -system {params.ini} \
      -threads {threads} \
      -somatic \
      --log {log}
    """