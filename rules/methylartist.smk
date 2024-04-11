rule methylartist:
    input:
        bam = "pacbio/alignment_{ref}/{sample}/{sample}.{ref}.haplotagged.bam",
        gtf = lambda wc: annot[wc.ref],
        ref = lambda wc: ref[wc.ref],
    output:
        directory("methylplot/pacbio_{ref}/{sample}/{roi}")
    conda:
        "env/methylartist.yml"
    params:
        roi = lambda wc: methylplot_roi[wc.roi][wc.ref],
        motif = "CG"
    log:
        "logs/{sample}_{ref}_methylartist_{roi}.log"
    shell:
        """
        mkdir -p {output}
        cd {output}
        methylartist locus \
            -i {params.roi} \
            --bams ../../../../{input.bam} \
            -g {input.gtf} \
            --ref {input.ref} \
            --phased \
            --color_by_hp \
            --ignore_ps \
            --motif {params.motif} \
            > ../../../../{log} 2>&1
        """