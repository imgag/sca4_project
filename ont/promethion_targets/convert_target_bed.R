library(data.table)

setwd("W:/users/ahgrosc1/data/promethion_targets")

dt <- fread("targeted-genes-prom.csv")
fai <- fread("W:/share/data/GRCh38/genomes/GRCh38.fa.fai", col.names = c("chr", "length", "offset", "linebases", "linewidth"))

print(dt)

dt[,c("chr", "coord") := tstrsplit(coordinates, ":"),]
dt[,c("start", "end") := tstrsplit(gsub(",", "", coord), "-"),]

dt2 <- dt[,.(chr, start = max(1, as.numeric(start) - 1e5), end = as.numeric(end) + 1e5, name),]


# Make sure that length is not higher then max chr length
e <- fai[dt2, on = "chr"][end>length,,]
dt2[name %in% e$name, end := e$length,]

fwrite(dt2, "prom_targets_220728.bed", sep = "\t", col.names = F)
