library(data.table)
library(ggplot2)

files <- Sys.glob("ont/Sample_*/repeat_expansions/straglr.tsv")
names <- fread("doc/names_ont.tsv")

dt_l <- list()

roi <- c(chrom = "chr16", start = 72787694, end = 72787756)

for (f in files) {
    sample = substr(f, 12, 22)
    dt <- fread(f)
    dt <- dt[`#chrom` == roi['chrom'] & start == roi['start'] & end == roi['end'] ,,]
    dt[, sample_gsvar := factor(sample, labels=names[id_folder == sample,id_pub,]),]
    dt_l[[sample]] <- dt
}

dt <- rbindlist(dt_l)

dt[allele<30, expanded := "no",]
dt[allele>=30, expanded := "yes"]

p <- ggplot(dt, aes(x = allele, y = copy_number, fill = expanded)) +
    #geom_point() +
    geom_boxplot() +
    facet_wrap(~sample_gsvar, nrow = 1, scales = "free_x") +
    theme_bw() +
    scale_fill_manual(values=c( "cornflowerblue", "coral3"))
p
ggsave("ont/ont_repeat_expansions.png", p,width = 18, height = 6)
ggsave("ont/ont_repeat_expansions.pdf", p,width = 18, height = 6)
