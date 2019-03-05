
library(biomaRt)
library(Biostrings)

animals <- list("human","chicken","pig","dog","mouse","rat")
animals_db <- list("hsapiens_gene_ensembl","ggallus_gene_ensembl","sscrofa_gene_ensembl","cfamiliaris_gene_ensembl","mmusculus_gene_ensembl","rnorvegicus_gene_ensembl")



#get all protein-coding gene IDs list for human

for(i in 1:6)
{
  animal <- animals[[i]]
  db <- animals_db[[i]]
  
  ensembl= useEnsembl(biomart="ensembl", dataset=db, version=88)
  
  proteinCodingGenes <- getBM(mart = ensembl, 
                              attributes = c("ensembl_gene_id","ensembl_transcript_id"),
                              filters = c("biotype"), 
                              values = list("protein_coding"))
  
  
  
  colnames(proteinCodingGenes) <- c("EnsembleGeneID","TranscriptID")
  
  write.csv(proteinCodingGenes, file = paste(animal,"_gene_transcript.csv"),row.names=FALSE)
  
  totaltranscriptsensemblv88<-unique(proteinCodingGenes$"TranscriptID")
  #write.table(totalgenesensemblv88, file= "totalgenesensemblv88.tab", sep = "\t", row.names = F)
  # extract gne nucleotide sequences
  
  gene <- getSequence(mart = ensembl, type = "ensembl_transcript_id", id = totaltranscriptsensemblv88, seqType = "cdna")
  
  
  colnames(gene) <- c("Sequence", "TranscriptGeneID")
  gene  <- gene[gene$Sequence != "Sequence unavailable", ]
  write.table(gene , file= "gene.tab", sep = "\t", row.names = F)
  
  #calculate the lenght for each sequence
  gene$length <- nchar(as.character(gene$Sequence))
  write.table(gene, file= "gene.tab", sep = "\t", row.names = F)
 
  
  write.csv(gene, file = paste(animal,"cdna_2.csv"),row.names=FALSE)
}

proteinCodingGenes <- getBM(mart = ensembl, 
                                  attributes = c("ensembl_gene_id","ensembl_transcript_id"),
                                  filters = c("biotype"), 
                                  values = list("protein_coding"))



colnames(proteinCodingGenes) <- c("EnsembleGeneID","TranscriptID")

write.csv(proteinCodingGenes, file = "MyData.csv",row.names=FALSE)

totaltranscriptsensemblv88<-unique(proteinCodingGenes$"TranscriptID")
#write.table(totalgenesensemblv88, file= "totalgenesensemblv88.tab", sep = "\t", row.names = F)
# extract gne nucleotide sequences

gene <- getSequence(mart = ensembl, type = "ensembl_transcript_id", id = totalgenesensemblv88, seqType = "cdna")


colnames(gene) <- c("Sequence", "EnsembleGeneID")
gene  <- gene[gene$Sequence != "Sequence unavailable", ]
write.table(gene , file= "gene.tab", sep = "\t", row.names = F)

#calculate the lenght for each sequence
gene$length <- nchar(as.character(gene$Sequence))
write.table(gene, file= "gene.tab", sep = "\t", row.names = F)

dna <- DNAStringSet(gene$Sequence)
names(dna) <- gene$"EnsembleGeneID"

write.csv(gene, file = "MyData.csv",row.names=FALSE)
# write to a file
#writeXStringSet(dna, filepath = "~/human.fas")

