suppressMessages(library(tidyverse))
suppressMessages(library(vcfR))

# Define default in/out files for testing if not run inside snakemake
if (exists("snakemake")) {
  input_vcf <- snakemake@input$vcf
#  input_vcf <- "pacbio/repeats_trgt/22018LRa001.repeat_catalog.GRCh38.vcf.gz"
#  out_tsv <- "pacbio/repeats_trgt/22018LRa001.repeat_catalog.expanded_loci.tsv"
#  out_folder <- "pacbio/repeats_trgt/22018LRa001.repeat_catalog.plot"
  out_tsv <- snakemake@output$tsv
  out_prefix <- str_replace(out_tsv, ".expanded_loci.tsv", "")
  out_folder <- snakemake@output$folder
} else {
  args = commandArgs(trailingOnly=TRUE)
  input_vcf <- args[1]
  out_prefix <- args[2]
}

plot_overview <- TRUE
create_folders <- TRUE

if(create_folders) dir.create(out_folder, recursive = TRUE)
print("STARTING SCRIPT")

#_____________________________________________________
# Read and format data from VCF File
vcf <- read.vcfR(input_vcf)
# Print INFO and FORMAT field descriptions and types
metaINFO2df(vcf, field = "INFO")
metaINFO2df(vcf, field = "FORMAT")

tb <- bind_cols(
    as_tibble(getFIX(vcf)),
    as_tibble(INFO2df(vcf)),
    as_tibble_col(extract.gt(vcf, element="GT")[,1], column_name = "GT"),
    as_tibble(str_split(extract.gt(vcf, element="AL"), ",", simplify = T)) %>% set_names(c("AL_1", "AL_2")),
    as_tibble(str_split(extract.gt(vcf, element="ALLR"), ",", simplify = T)) %>% set_names(c("ALLR_1", "ALLR_2")),
    as_tibble(str_split(extract.gt(vcf, element="SD"), ",", simplify = T)) %>% set_names(c("SD_1", "SD_2")),
    as_tibble(str_split(extract.gt(vcf, element="MC"), ",", simplify = T)) %>% set_names(c("MC_1", "MC_2")),
    as_tibble(str_split(extract.gt(vcf, element="MS"), ",", simplify = T)) %>% set_names(c("MS_1", "MS_2")),
    as_tibble(str_split(extract.gt(vcf, element="AP"), ",", simplify = T)) %>% set_names(c("AP_1", "AP_2")),
    as_tibble(str_split(extract.gt(vcf, element="AM"), ",", simplify = T)) %>% set_names(c("AM_1", "AM_2")),
) %>%
    mutate_at('POS', as.integer) %>%
    mutate_at('QUAL', as.integer) %>%
    mutate_at('AL_1', as.integer) %>%
    mutate_at('AL_2', as.integer) %>%
    mutate_at('SD_1', as.integer) %>%
    mutate_at('SD_2', as.integer) %>%
    mutate_at('MC_1', as.integer) %>%
    mutate_at('MC_2', as.integer) %>%
    mutate_at('AM_1', as.integer) %>%
    mutate_at('AM_2', as.integer) %>%
    mutate('SAMPLE' = colnames(vcf@gt)[-1]) %>%
    mutate(REF_LENGTH = nchar(REF))

#_____________________________________________________
# Identify expanded sites

MIN_EXPANSION_BP <- 50              # ALT must be at least 50bp longer then REF
MIN_EXPANSION_FACTOR <- 3           # or ALT must be at least n times the REF length
MIN_SUPPORTING_READS <- 3           # MIN number of reads per allele

tvf <- tb %>%
    filter(
        AL_1 > MIN_EXPANSION_BP + REF_LENGTH     |
        AL_2 > MIN_EXPANSION_BP + REF_LENGTH     |
        AL_1 > MIN_EXPANSION_FACTOR * REF_LENGTH |
        AL_2 > MIN_EXPANSION_FACTOR * REF_LENGTH &
        SD_1 >= MIN_SUPPORTING_READS             &
        SD_2 >= MIN_SUPPORTING_READS
    )

print(paste("Found", nrow(tvf), "expanded sites"))

write_delim(tvf, paste0(out_prefix, ".expanded_loci.tsv"), delim="\t")

#_____________________________________________________
# Write folders for expanded sites

if (create_folders) {
    for (i in seq_along(tvf$TRID)) {
        row <- tvf[i,]
        outfile <- paste0(out_prefix, ".plot/", row$TRID,"/repeat_locus_info.tsv")
        dir.create(file.path(out_folder, row$TRID), recursive = TRUE)
        write_delim(row, file = outfile, delim = "\t")
    }
}

#____________________________________________________

MAX_LOCI <- 400
MAX_NAMED_LOCI <- 60

# Reduce set (favoring EXPANDED) for plotting.
tb2 <- tb %>% 
    mutate(EXPANDED = if_else(TRID %in% tvf$TRID, "Expanded", "Normal")) %>%
    arrange(EXPANDED, desc((AL_1 + AL_2) / 2* REF_LENGTH)) %>%
    head(MAX_LOCI)

# Long table format for plotting
tb3 <- tb2 %>%
    pivot_longer(
        c(ends_with("_1"), ends_with("_2")),
        names_pattern = "(.*)_(.)",
        names_to = c(".value", "ALLELE")
    )%>%
    mutate(
        ALLR_lower = as.integer(str_split(ALLR, "-", simplify =T)[,1]),
        ALLR_upper = as.integer(str_split(ALLR, "-", simplify =T)[,2]),
    )

p1 <- ggplot(tb3, aes(x = TRID, colour = ALLELE, shape = EXPANDED)) +
    geom_pointrange(aes(y = AL, ymin = ALLR_lower, ymax = ALLR_upper)) +
    geom_point(aes(y = REF_LENGTH), colour = "grey40") +
    theme_classic() +
    scale_shape_manual(values = c(8,19)) +
    labs(
        x = "Repeat loci",
        y = "Alternative allele length" ) 

if (n_distinct(tb3$TRID) > MAX_NAMED_LOCI) {
    p1 <- p1 + theme(axis.text.x = element_blank())
} else {
    p1 <- p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

ggsave(paste0(out_prefix, ".repeat_loci_length.png"), plot = p1, width = 14, height = 6)
ggsave(paste0(out_prefix, ".repeat_loci_length.pdf"), plot = p1, width = 14, height = 6)

p2 <- ggplot(tb3, aes(x = TRID, col = ALLELE, shape = EXPANDED)) +
    geom_pointrange(aes(y = AL/REF_LENGTH, ymin = ALLR_lower/REF_LENGTH, ymax = ALLR_upper/REF_LENGTH)) +
    geom_point(aes(y = 1), col = "grey40") + 
    scale_y_log10() +
    theme_classic() +
    scale_shape_manual(values = c(8,19)) +
    labs(
        x = "Repeat loci",
        y = "Alternative allele length / Ref length" ) 

if (n_distinct(tb3$TRID) > MAX_NAMED_LOCI) {
    p2 <- p2 + theme(axis.text.x = element_blank())
} else {
    p2 <- p2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

ggsave(paste0(out_prefix, ".repeat_loci_normalized.png"), plot = p2, width = 14, height = 6)
ggsave(paste0(out_prefix, ".repeat_loci_normalized.pdf"), plot = p2, width = 14, height = 6)