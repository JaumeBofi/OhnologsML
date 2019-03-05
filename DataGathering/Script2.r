####################################################################
#human data source
ensembl= useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=88)

#get all protein-coding gene IDs list for human

proteinCodingTranscripts <- getBM(mart = ensembl, 
                                  attributes = c("ensembl_gene_id"),
                                  filters = c("biotype"), 
                                  values = list("protein_coding"))
colnames(proteinCodingGenes) <- c("EnsembleGeneID")

totalgenesensemblv88<-unique(proteinCodingGenes$EnsembleGeneID)
#write.table(totalgenesensemblv88, file= "totalgenesensemblv88.tab", sep = "\t", row.names = F)
# extract gne nucleotide sequences
gene <- getSequence(mart = ensembl, type = "ensembl_gene_id", id = totalgenesensemblv88, seqType = "gene_exon_intron")
colnames(gene) <- c("Sequence", "EnsembleGeneID")
gene  <- gene[gene$Sequence != "Sequence unavailable", ]
write.table(gene , file= "gene.tab", sep = "\t", row.names = F)

#calculate the lenght for each sequence
gene$length <- nchar(as.character(gene$Sequence))
write.table(gene, file= "gene.tab", sep = "\t", row.names = F)

dna <- DNAStringSet(gene$Sequence)
names(dna) <- gene$EnsembleGeneID

# write to a file
writeXStringSet(dna, filepath = "~/jaume/geneNucleoSeq_hsv88.fas")