############# 

proteinCodingGenes <- getBM(mart = ensembl,

                                  attributes = c("ensembl_gene_id"),

                                  filters = c("biotype"),

                                  values = list("protein_coding"))

colnames(proteinCodingGenes) <- c("EnsembleGeneID")

 

totalgenesensemblv88<-unique(proteinCodingGenes$EnsembleGeneID)

#write.table(totalgenesensemblv88, file= "totalgenesensemblv88.tab", sep = "\t", row.names = F)


gene <- getSequence(mart = ensembl, type = "ensembl_gene_id", id = totalgenesensemblv88, seqType = "gene_exon_intron")

colnames(gene) <- c("Sequence", "EnsembleGeneID")

gene  <- gene[gene$Sequence != "Sequence unavailable", ] # esta linea es para que solo extraiga las secuencias de genes disponibles en la base de dato de ensembl 

write.table(gene , file= "gene.tab", sep = "\t", row.names = F)

#################

