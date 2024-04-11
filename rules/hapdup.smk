#_______ SV HAPDUP DIPDIFF _________________________________________________________#

rule map_to_assembly:
    input:
        ref = "pacbio/assembly_hifiasm/{sample}/{sample}.asm_hifi.bp.p_ctg.fasta",
        xml = "pacbio/{sample}.ccs.xml"
    output:
        bam = "pacbio/alignment_assembly/{sample}/{sample}.asm_hifi.bam"
    conda:
        "env/pbmm2.yml"
    log:
        "logs/{sample}_pbmm2.log"
    threads:
        16
    shell:
        """
        pbmm2 align --preset CCS --unmapped -j {threads} {input.ref} {input.xml} {output} --sort -j4 -J 2 >{log} 2>&1
        """

rule hapdup:
  input:
    bam = "pacbio/alignment_assembly/{sample}/{sample}.asm_hifi.bam",
    fa = "pacbio/assembly_hifiasm/{sample}/{sample}.asm_hifi.bp.p_ctg.fasta"
  output:
    f1 = "pacbio/phased_assembly/{sample}_hapdup/hapdup_dual_1.fasta",
    f2 = "pacbio/phased_assembly/{sample}_hapdup/hapdup_dual_2.fasta"
  log:
    "logs/{sample}_hapdup.log"
  threads:
    24
  resources:
    queue="research_srv023"
  params:
    image = "caspargross/hapdup:latest" 
#    singularity exec \
#      --bind "$(realpath $(dirname {input.bam}))":"/mnt/alignment" \
#      --bind "$(realpath $(dirname {input.fa}))":"/mnt/assembly" \
#      --bind "$(realpath $(dirname {output.f1}))":"/mnt/hapdup" \
  shell:
    """
    docker run \
      -v "$(realpath $(dirname {input.bam}))":"/mnt/alignment" \
      -v "$(realpath $(dirname {input.fa}))":"/mnt/assembly" \
      -v "$(realpath $(dirname {output.f1}))":"/mnt/hapdup" \
      -u `id -u`:`id -g` \
      {params.image} \
      hapdup \
      --assembly /mnt/assembly/$(basename {input.fa})\
      --bam /mnt/alignment/$(basename {input.bam}) \
      --out-dir /mnt/hapdup \
      -t {threads} \
      --rtype hifi \
      > {log} 2>&1
   """

rule hapdiff:
  input:
    ref = lambda w: ref[w.ref],
    f1 = "pacbio/phased_assembly/{sample}_hapdup/hapdup_dual_1.fasta",
    f2 = "pacbio/phased_assembly/{sample}_hapdup/hapdup_dual_2.fasta"
  output:
    vcf = "pacbio/variant_calling/{sample}_{ref}_dipdiff/hapdiff_phased.vcf.gz",
    vcfu = "pacbio/variant_calling/{sample}_{ref}_dipdiff/hapdiff_unphased.vcf.gz"
  log:
    "logs/{sample}_{ref}_dipdiff.log"
  threads:
    12
  resources:
    queue="research_srv023"
  params:
#    image = "/mnt/users/ahgrosc1/container/hapdiff_0.7.sif"
    image = "mkolmogo/hapdiff:0.7"
  shell: 
#    singularity exec \
#        --bind "$(realpath $(dirname {input.ref}))":"/mnt/ref" \
#        --bind "$(realpath $(dirname {input.f1}))":"/mnt/hapdup" \
#        --bind "$(realpath $(dirname {output.vcf}))":"/mnt/output" \
    """
      docker run \
        -v "$(realpath $(dirname {input.ref}))":"/mnt/ref" \
        -v "$(realpath $(dirname {input.f1}))":"/mnt/hapdup" \
        -v "$(realpath $(dirname {output.vcf}))":"/mnt/output" \
        -u `id -u`:`id -g` \
        {params.image}\
        hapdiff.py \
        --reference "/mnt/ref/$(basename {input.ref})" \
        --pat "/mnt/hapdup/hapdup_dual_1.fasta" \
        --mat "/mnt/hapdup/hapdup_dual_2.fasta" \
        --out-dir /mnt/output \
        -t {threads} \
        >{log} 2>&1
    """

