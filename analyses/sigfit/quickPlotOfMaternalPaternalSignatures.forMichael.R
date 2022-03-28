# quick plot of different signatures for kh and michael
plotdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/sigfit/comparisonOfAgingSignaturesForMichael/"
dir.create(plotdir)
aging_sigs_manyVersions <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/information/jonsson_parental_aging/jonsson.agingsignatures.tables9.convertedtoproportionsinexcel.REVCOMP.REORDEREDTOMATCHMYDATA.LOTSOFDIFFERENTVERSIONSFORPLOTTING.dontuseforfitting.txt",header=T)
head(aging_sigs_manyVersions)

aging_sigs_manyVersions_melt <- melt(aging_sigs_manyVersions)
head(aging_sigs_manyVersions_melt)
ggplot(aging_sigs_manyVersions_melt,aes(x=variable,y=value,fill=signature))+
  geom_col(position="dodge")+
  scale_fill_brewer(palette = "Set1")

require(sigfit)
# need rownames but no signature label:
rownames(aging_sigs_manyVersions) <- aging_sigs_manyVersions$signature
# get rid of signature column
aging_sigs_manyVersions <- select(aging_sigs_manyVersions,-"signature")

GetCosineSimilarity_betweenSignatures <- function(de_novo_signature_set,a_priori_signature_set){
  # make sure names are in same order
  if(sum(colnames(de_novo_signature_set)!=colnames(a_priori_signature_set))!=0){
    print("different column names for signature sets")
    break
  }
  nsig_denovo=dim(de_novo_signature_set)[1] # number of sigs
  nsig_apriori=dim(a_priori_signature_set)[1] # number of sigs
  allCos=data.frame()
  for(i in seq(1,nsig_denovo)){
    for(j in seq(1,nsig_apriori)){
      cossimilarity=cosine_sim(as.numeric(de_novo_signature_set[i,]),as.numeric(a_priori_signature_set[j,])) # get cosine similarity (this is how sigfit does it in plotting code -- not very efficient, but oh well )
      cosdf = data.frame(de_novo_signature=rownames(de_novo_signature_set)[i],a_priori_signature=rownames(a_priori_signature_set)[j],cosine_similarity=cossimilarity)
      allCos=bind_rows(allCos,cosdf)
    }}
  return(allCos)
}


cosineSimilaritiesALlVersionsOfAgingSIgnatures <- GetCosineSimilarity_betweenSignatures(aging_sigs_manyVersions,aging_sigs_manyVersions)

p1 <-ggplot(cosineSimilaritiesALlVersionsOfAgingSIgnatures,aes(x=de_novo_signature,y=a_priori_signature,fill=cosine_similarity))+
  geom_tile()+
  scale_fill_viridis_c()+
  geom_text(aes(label=round(cosine_similarity,2)))+
  ylab("")+
  xlab("")+
  ggtitle('cosine similarity')+
  theme(axis.text.x=element_text(angle=45,hjust=1))

p1
ggsave(paste0(plotdir,"cosinesimplot.png"),p1,height=8,width=14)
# get clr-aitchison distance between signatures :: try this!
GetAitchisonDistance_betweenSignatures <- function(signature_set_1,signature_set_2){
  # make sure names are in same order
  if(sum(colnames(signature_set_1)!=colnames(signature_set_2))!=0){
    print("different column names for signature sets")
    break
  }
  nsig1=dim(signature_set_1)[1] # number of sigs
  nsig2=dim(signature_set_2)[1] # number of sigs
  allDist=data.frame()
  for(i in seq(1,nsig1)){
    for(j in seq(1,nsig2)){
      # need to make a matrix: 
      dist=tidy(dist(rbind(clr(as.numeric(signature_set_1[i,])),clr(as.numeric(signature_set_2[j,]))))) # get cosine similarity (this is how sigfit does it in plotting code -- not very efficient, but oh well )
      distdf = data.frame(signature_set_1=rownames(signature_set_1)[i],signature_set_2=rownames(signature_set_2)[j],aitchison_clr_dist=dist$distance)
      allDist=bind_rows(allDist,distdf)
    }}
  return(allDist)
  
  
}

aitchisonDistAllVersionsOfAgingSIgnatures <- GetAitchisonDistance_betweenSignatures(aging_sigs_manyVersions,aging_sigs_manyVersions)


p2 <- ggplot(aitchisonDistAllVersionsOfAgingSIgnatures,aes(x=signature_set_1,y=signature_set_2,fill=aitchison_clr_dist))+
  geom_tile()+
  scale_fill_viridis_c(direction=-1)+
  geom_text(aes(label=round(aitchison_clr_dist,2)))+
  ylab("")+
  xlab("")+
  ggtitle("aitchison distance")+
  theme(axis.text.x=element_text(angle=45,hjust=1))
p2
ggsave(paste0(plotdir,"aitchisondist.png"),p2,height=8,width=14)

