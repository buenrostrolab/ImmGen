homer <- read_tsv("../data/homer_short.bed", col_names = FALSE)
homer <- data.frame(homer)
homer$match <- paste0(homer[,1], ":", homer[,3])

real <- read.table("../data/161021_peaks_quick.bed")
real$match <- paste0(real[,1], ":", real[,3])
real$id  <- 1:nrow(real)
merged <- merge(real, homer, by.x = "match", by.y = "match")
merged <- merged[order(merged$id), ]
m <- merged[,c(2,3,4,10,11,13,5)]


rownames(m) <- NULL
colnames(m) <- NULL

m[,4] <- grepl("promoter", m[,4])
write.table(m, file = "../data/161021_peaks_annotated.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)



# Done
write.table(data.frame(data.matrix(motif_score)), file = "../data/161021_peaksXscores.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Need to do
write.table(data.frame(data.matrix(motif_ix)), file = "../data/161021_peaksXmatch.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(data.matrix((deviations$z))), file = "../data/161021_motifsXsamples_zscores.tsv", sep = "\t", row.names = TRUE, quote = FALSE)

