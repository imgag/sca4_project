rule fix_gvcf_samplename:
    input:
        gvcf = "pacbio/variant_calling/{sample}_{ref}_deepvariant/{sample}.dv.gvcf.gz"
    output:
        gvcf = "pacbio/variant_calling/{sample}_{ref}_deepvariant/{sample}.fixed.dv.gvcf.gz"
    conda:
        "env/glnexus.yml"
    shell:
        """
        bcftools reheader -s <(echo {wildcards.sample}) -o {output.gvcf} {input.gvcf}
        """

rule merge_snv_gvcf:
    input:
        gvcf = expand("pacbio/variant_calling/{sample}_{{ref}}_deepvariant/{sample}.fixed.dv.gvcf.gz", sample = list(samples_pb.keys()))
    output:
        gvcf = "pacbio/merged_variants/all_snv-indel_dv.{ref}.gvcf"
    conda:
        "env/glnexus.yml"
    log:
        "logs/merge_snv_glnexus_{ref}.log"
    threads:
        6
    shell:
        """
        glnexus_cli \
        --config DeepVariantWGS \
        --threads {threads} \
        --mem-gbytes 200G \
        --dir pacbio/merged_variants/{wildcards.ref}.GLnexus.DB \
        {input.gvcf} > {output.gvcf} 2>{log} 
        """

rule zip_merge_snv_gvcf:
    input:
        gvcf = "pacbio/merged_variants/all_snv-indel_dv.{ref}.gvcf"
    output:
        gvcf = "pacbio/merged_variants/all_snv-indel_dv.{ref}.gvcf.gz"
    conda:
        "env/glnexus.yml"
    log:
        "logs/bgzip_merge_snv_glnexus_{ref}.log"
    shell:
        """
        bcftools view {input.gvcf} \
        | bgzip -c > \
        {output.gvcf}
        tabix {output.gvcf}
        """

rule unzip_survivor_input:
    input:
        vcf = "pacbio/variant_calling/{sample}_{ref}_dipdiff/hapdiff_unphased.vcf.gz"
    output:
        vcf = "pacbio/variant_calling/{sample}_{ref}_dipdiff/hapdiff_unphased.vcf"
    shell:
        """
        gunzip -c {input} > {output}
        """

rule merge_sv_survivor:
    input:
        vcf = expand("pacbio/variant_calling/{sample}_{{ref}}_dipdiff/hapdiff_unphased.vcf", sample = samples_pb_reduced)
    output:
        vcf = "pacbio/merged_variants/all_sv_assembly.{ref}.vcf"
    conda:
        "env/survivor.yml"
    log:
        "logs/merge_sv_{ref}_survivor.log"
    # SURVIVOR options
    params:
        d = 1000,   # maximum allowed distance between SVs
        c = 1,      # minimum number of supporting callers
        t = 1,      # variants need to agree on type (=1),
        s = 1,      # variants need to agree on strand (=1),
        x = 0,      # disabled parameter
        l = 50,     # minimum length of variants
    shell:
        """
        SURVIVOR merge \
        <(echo {input.vcf} | tr -s '[:blank:]' '\n') \
        {params.d} {params.c} {params.t} {params.s} {params.x} {params.l} \
        {output.vcf}  2> {log}
        """