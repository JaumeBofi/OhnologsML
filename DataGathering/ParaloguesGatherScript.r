
library(biomaRt)
ensembl <- useMart("ensembl")

filepath = "E:\StuffRemoteStorage\Thesis\Ohnologs\\"
filepath_transcripts <- paste(filepath,animal,"\\",sep="")
filepath_ohnologs <- paste(filepath_transcripts,level,"\\",sep="")
filepath_ohnologs


ohnolohgs <- scan(paste(filepath_ohnologs,"HUMAN.Ohnologs-List-Complete.txt",sep=""), what="", sep=",")


ensembl_human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

attributes = listAttributes(ensembl_human)

getBM(attributes = c('ensembl_gene_id','hsapiens_paralog_ensembl_gene'),
      filters = 'ensembl_gene_id',
      values = c(ohnolohgs),
      mart = ensembl_human)





library(biomaRt)
ensembl <- useMart("ensembl")

animalLists <- list("Human","Chicken","Pig","Dog","Mouse","Rat")
levelLists <- list("Strict","Relaxed","Intermediate")
animals_db <- list("hsapiens_gene_ensembl","ggallus_gene_ensembl","sscrofa_gene_ensembl","cfamiliaris_gene_ensembl","mmusculus_gene_ensembl","rnorvegicus_gene_ensembl")
animals_family <- list("hsapiens","ggallus","sscrofa","cfamiliaris","mmusculus","rnorvegicus")
filepath = "E:\\StuffRemoteStorage\\Thesis\\Ohnologs\\"
paralogues <- NULL

for(i in 1:6)
{
  animal <- animalLists[[i]]
  animal_db <- animals_db[[i]]
  animal_family <- animals_family[[i]]
  
  print(paste("Animal-",animal))
  for(j in 1:3)
  {
    level <- levelLists[[j]]
    print(paste("Level-",level))  
    
    filepath_transcripts <- paste(filepath,animal,"\\",sep="")
    filepath_ohnologs <- paste(filepath_transcripts,level,"\\",sep="")
    
    
    ohnolohgs <- scan(paste(filepath_ohnologs,paste(toupper(animal),".Ohnologs-List-Complete.txt",sep=""),sep=""), what="", sep=",")
    ensembl_human <- useMart("ensembl", dataset = animal_db)
    paralogues <- getBM(attributes = c('ensembl_gene_id',paste(animal_family,'_paralog_ensembl_gene',sep="")),
          filters = 'ensembl_gene_id',
          values = c(ohnolohgs),
          mart = ensembl_human)
    
    write.csv(paralogues, file = paste(filepath_ohnologs,paste(toupper(animal),".Paralogues-Non-Filter.csv",sep=""),sep=""),row.names=FALSE)
  }
}