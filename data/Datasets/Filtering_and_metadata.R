library(dplyr)
`%nin%` = Negate(`%in%`)

listado_variantes.tab <-read.table("listado_variantes.txt")
archivos <- listado_variantes.tab$V2 
listado_variantes.tab$V2 <- gsub(".gbk","",archivos)

metadatos_completos <- read.table(file = 'metadata.tsv', sep = '\t', header = TRUE, na.strings = "", fill=TRUE,quote = "")

metadatos_muestra <- merge(metadatos_completos,listado_variantes.tab, by.x = 2, by.y= 2)
colnames(metadatos_muestra)[18]<-"Variant"

#Variantes
variants_of <- c("alpha","beta","gamma","delta","eta","iota","kappa","lambda")

#Verificando linajes
#alpha
alpha.tab <- metadatos_muestra[which(metadatos_muestra$Variant=="alpha"),]
nosonalpha <- alpha.tab[alpha.tab$Lineage != "B.1.1.7",]
remover_alpha <- nosonalpha$Accession.ID

corrected_alpha.tab <- alpha.tab[alpha.tab$Lineage == "B.1.1.7",]
corrected_alpha.tab <- na.omit(corrected_alpha.tab)
rowstoremove <- c()
for (i in remover_alpha){
  system(paste("rm alphas/*",i,".gbk",sep = ""))
  rowtoremove <- which(alpha.tab$Accession.ID==i)
  rowstoremove<-append(rowstoremove, rowtoremove, after = length(rowstoremove))
}

alpha.tab <- alpha.tab[-rowstoremove,]

#beta
beta.tab <- metadatos_muestra[which(metadatos_muestra$Variant=="beta"),]
nosonbeta <- beta.tab[beta.tab$Lineage != "B.1.351",]

#delta
delta.tab <- metadatos_muestra[which(metadatos_muestra$Variant=="delta"),]
nosondelta <- delta.tab[delta.tab$Lineage != "B.1.617.2",]

#epsilon
epsilon.tab <- metadatos_muestra[which(metadatos_muestra$Variant=="epsilon"),]
nosonepsilon <- epsilon.tab[epsilon.tab$Lineage != "B.1.427",]

#eta
eta.tab <- metadatos_muestra[which(metadatos_muestra$Variant=="eta"),]
nosoneta <- eta.tab[eta.tab$Lineage != "B.1.525",]

#gamma
gamma.tab <- metadatos_muestra[which(metadatos_muestra$Variant=="gamma"),]
nosongamma <- gamma.tab[gamma.tab$Lineage != "P.1",]

#iota
iota.tab <- metadatos_muestra[which(metadatos_muestra$Variant=="iota"),]
nosoniota <- iota.tab[iota.tab$Lineage != "B.1.526",]

#kappa
kappa.tab <- metadatos_muestra[which(metadatos_muestra$Variant=="kappa"),]
nosonkappa <- kappa.tab[kappa.tab$Lineage != "B.1.617.1",]

#lamda
lambda.tab <- metadatos_muestra[which(metadatos_muestra$Variant=="lambda"),]
nosonlambda <- lambda.tab[lambda.tab$Lineage != "C.37",]

#zeta
zeta.tab <- metadatos_muestra[which(metadatos_muestra$Variant=="zeta"),]
nosonzeta <- zeta.tab[zeta.tab$Lineage != "P.2",]

#BUILDING METADATA FILES
#VOC
VOCdelete <- c("FR989824","FR989826","HG999984")
VOC_metadata <- rbind(alpha.tab,beta.tab,gamma.tab,delta.tab)
VOC_metadata <- VOC_metadata[which(VOC_metadata$Accession.ID %nin% VOCdelete ),]
for (i in unique(VOC_metadata$Variant)){
  print(paste((length(which(VOC_metadata$Variant==i))),i))
}

#VOI
VOI_metadata <- rbind(eta.tab,iota.tab,kappa.tab,lambda.tab)
for (i in unique(VOI_metadata$Variant)){
  print(paste((length(which(VOI_metadata$Variant==i))),i))
}

#OTHERS
others_listado <- read.table("lista_others.txt")
Others_metadata <- merge(others_listado,metadatos_muestra,by.x=1,by.y = 1)

#WRITE METADATA FILES
write.table(VOC_metadata, file='VOC_metadata.tsv', quote=FALSE, sep='\t', col.names = NA)
write.table(VOI_metadata, file='VOI_metadata.tsv', quote=FALSE, sep='\t', col.names = NA)
write.table(Others_metadata, file='Others_metadata.tsv', quote=FALSE, sep='\t', col.names = NA)

#REMOVE UNANNOTATED FILES
Folders <- c("alphas","betas","gammas","deltas","etas","iotas","kappas","lambdas","others")

#alphas
alphas_translation_grep <- as.data.frame(system('grep -c translation alphas/*.gbk',intern=TRUE))
for (i in 1:length(alphas_translation_grep[,1])){
  alphas_translation_grep$is_empty[i]<- grepl(":0",alphas_translation_grep[i,1] , fixed = TRUE)
}
#betas
betas_translation_grep <- as.data.frame(system('grep -c translation betas/*.gbk',intern=TRUE))
for (i in 1:length(betas_translation_grep[,1])){
  betas_translation_grep$is_empty[i]<- grepl(":0",betas_translation_grep[i,1] , fixed = TRUE)
}
#gammas
gammas_translation_grep <- as.data.frame(system('grep -c translation gammas/*.gbk',intern=TRUE))
for (i in 1:length(gammas_translation_grep[,1])){
  gammas_translation_grep$is_empty[i]<- grepl(":0",gammas_translation_grep[i,1] , fixed = TRUE)
}
#deltas
deltas_translation_grep <- as.data.frame(system('grep -c translation deltas/*.gbk',intern=TRUE))
for (i in 1:length(deltas_translation_grep[,1])){
  deltas_translation_grep$is_empty[i]<- grepl(":0",deltas_translation_grep[i,1] , fixed = TRUE)
}
#etas
etas_translation_grep <- as.data.frame(system('grep -c translation etas/*.gbk',intern=TRUE))
for (i in 1:length(etas_translation_grep[,1])){
  etas_translation_grep$is_empty[i]<- grepl(":0",etas_translation_grep[i,1] , fixed = TRUE)
}
#iotas
iotas_translation_grep <- as.data.frame(system('grep -c translation iotas/*.gbk',intern=TRUE))
for (i in 1:length(iotas_translation_grep[,1])){
  iotas_translation_grep$is_empty[i]<- grepl(":0",iotas_translation_grep[i,1] , fixed = TRUE)
}
#kappas
kappas_translation_grep <- as.data.frame(system('grep -c translation kappas/*.gbk',intern=TRUE))
for (i in 1:length(kappas_translation_grep[,1])){
  kappas_translation_grep$is_empty[i]<- grepl(":0",kappas_translation_grep[i,1] , fixed = TRUE)
}
#lambdas
lambdas_translation_grep <- as.data.frame(system('grep -c translation lambdas/*.gbk',intern=TRUE))
for (i in 1:length(lambdas_translation_grep[,1])){
  lambdas_translation_grep$is_empty[i]<- grepl(":0",lambdas_translation_grep[i,1] , fixed = TRUE)
}
#others
others_translation_grep <- as.data.frame(system('grep -c translation others/*.gbk',intern=TRUE))
for (i in 1:length(others_translation_grep[,1])){
  others_translation_grep$is_empty[i]<- grepl(":0",others_translation_grep[i,1] , fixed = TRUE)
}

#DELETE EM
colnames(alphas_translation_grep)<-c("file","is_empty")
colnames(betas_translation_grep)<-c("file","is_empty")
colnames(gammas_translation_grep)<-c("file","is_empty")
colnames(deltas_translation_grep)<-c("file","is_empty")
colnames(etas_translation_grep)<-c("file","is_empty")
colnames(iotas_translation_grep)<-c("file","is_empty")
colnames(kappas_translation_grep)<-c("file","is_empty")
colnames(lambdas_translation_grep)<-c("file","is_empty")
colnames(others_translation_grep)<-c("file","is_empty")
Files.db <- rbind(alphas_translation_grep,betas_translation_grep,gammas_translation_grep,deltas_translation_grep,etas_translation_grep,iotas_translation_grep,kappas_translation_grep,lambdas_translation_grep,others_translation_grep)
Empty_files.db <- Files.db[which(Files.db$is_empty == TRUE),]
Empty_files.db$filedel <- paste("rm",gsub(":0","",Empty_files.db$file),sep = " ")
Empty_files.db$Identifiers <- gsub(".gbk:0","",gsub("^.*_","",Empty_files.db$file))
#write.table(Empty_files.db, file='Empty_files.csv', quote=FALSE, sep=',', col.names = NA)
EMPTYFILES <- read.table(file = 'Empty_files.csv', sep = ',', header = TRUE, na.strings = "", fill=TRUE,quote = "")

#Delete rows of empty files from metadata
Filtered_VOC_metadata <- VOC_metadata[which(VOC_metadata$Accession.ID %nin% EMPTYFILES$Identifiers),]
Filtered_VOI_metadata <- VOI_metadata[which(VOI_metadata$Accession.ID %nin% EMPTYFILES$Identifiers),]
Filtered_Others_metadata <- Others_metadata[which(Others_metadata$Accession.ID %nin% EMPTYFILES$Identifiers),]

#WRITE FILTERED METADATA FILES
write.table(Filtered_VOC_metadata, file='VOC_metadata.tsv', quote=FALSE, sep='\t', col.names = NA)
write.table(Filtered_VOI_metadata, file='VOI_metadata.tsv', quote=FALSE, sep='\t', col.names = NA)
write.table(Others_metadata, file='Others_metadata.tsv', quote=FALSE, sep='\t', col.names = NA)






