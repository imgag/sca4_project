find ../ -name *repeat_catalog.GRCh38.vcf.gz -exec sh -c "zcat {} | grep 72787694 " \; > zfhx3_repeat_lengts.tsv
