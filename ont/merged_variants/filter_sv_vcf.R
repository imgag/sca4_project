library(tidyverse)
library(vcfR)
library(eulerr)

for (tech in c("dipdiff", "cutesv", "sniffles")){

    vcf <- read.vcfR(paste0("merged_sv.", tech, ".vcf"))


    #___________________________________________________
    # EXAMPLE COMMANDS

    # VCF overview
    #print(vcf)

    # Metadata information
    #queryMETA(vcf)

    # Get available fixed regions
    #head(getFIX(vcf),1)

    # GT region (per sample information)
   # head(vcf@gt[1:6, 1:4])

    # Extract information from the GT fields
    #head(extract.gt(vcf, element="GT"))


    #___________________________________________________
    # FORMAT DATA 

    samplenames <- colnames(vcf@gt)[-1]
    # Overwrite in this case
    #samplenames <- c("22018LRa001", "22018LRa002", "22018LRa003", "22018LRa004")
    #samplenames <- c("F1 father (affected)", "F1 mother", "F1 daughter", "F1 daughter (affected)")

    # Overall number of variants with supporting variants
    sup_vec <- extract.info(vcf, "SUPP_VEC")

    # Convert dataframe to tidy
    id <- (vcf@fix)[,'ID']                          # Need to make ID field unique by appending row index
    id_unique <- paste0(id, ".", 1:length(id))
    vcf@fix[,'ID'] <- id_unique

    tb <- vcfR2tidy(vcf)                            # Use vcfR function to convert to tibble
    head(tb)


    #___________________________________________________
    # PLOTs FOR WHOLE GENOME

    # Number ob supporting samples per SV types
    df <- INFO2df(vcf)

    ggplot(df, aes(x=SUPP, fill=SVTYPE)) +
        geom_histogram(stat="count", position = "dodge", colour="grey20") +
        theme_classic() +
        xlab("number of samples with variant")

    ggsave(paste0("plot/", tech, ".sv_support_freq.whole_genome.png"), width = 5, height = 4)

    # Plot set of supporting samples as euler diagram
    sup_mat <- t(sapply(sup_vec, 
        function(x) unlist(strsplit(x, "")) == "1",
        USE.NAMES=FALSE))
    colnames(sup_mat) <- samplenames

    fit_e <- euler(sup_mat)
    png(paste0("plot/", tech, ".sv_set.whole_genome.png"), width = 600, height = 600)
    plot(fit_e,
        quantities = TRUE,
        labels = list(font =4),
        lty = 1:4)
    dev.off()
    pdf(paste0("plot/", tech, ".sv_set.whole_genome.pdf"), width = 600, height = 600)
    plot(fit_e,
        quantities = TRUE,
        labels = list(font =4),
        lty = 1:4)
    dev.off()

    #___________________________________________________
    # PLOTs FOR TARGET REGION 

    # Number ob supporting samples per SV types
    tb <- bind_cols(
        as_tibble(getFIX(vcf)),
        as_tibble(INFO2df(vcf)),
        ) %>%
        mutate_at('POS', as.integer)

    tb_roi <- tb %>%
        filter(CHROM == "chr16") %>%
        filter(POS > 67788110) %>%
        filter(POS < 79817700)

    ggplot(tb_roi, aes(x=SUPP, fill=SVTYPE)) +
        geom_histogram(stat="count", position = "dodge", colour="grey20") +
        theme_classic() +
        xlab("number of samples with variant")

    ggsave(paste0("plot/", tech, ".sv_support_freq.region.png"), width = 5, height = 4)

    # Plot set of supporting samples as euler diagram
    sup_mat <- t(sapply(tb_roi$SUPP_VEC, 
        function(x) unlist(strsplit(x, "")) == "1",
        USE.NAMES=FALSE))
    colnames(sup_mat) <- samplenames

    fit_e_r <- euler(sup_mat)
    png(paste0("plot/", tech, ".sv_set.region.png"), width = 600, height = 600)
    plot(fit_e_r,
        quantities = TRUE,
        labels = list(font =4),
        lty = 1:4)
    dev.off()
    pdf(paste0("plot/", tech, ".sv_set.region.pdf"), width = 600, height = 600)
    plot(fit_e_r,
        quantities = TRUE,
        labels = list(font =4),
        lty = 1:4)
    dev.off()


    #___________________________________________________
    # EXPORT CANDIDATES

    tb_roi %>% filter(SUPP_VEC == "10011111") %>% write_delim(paste0(tech,".SV_candidates_region_F1.tsv"))
}
