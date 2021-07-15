
#LOADING LIBRARY
library("tidyverse")
library("ggrepel")

#MAKING OUTDIR
outdir <- "../results/Complexity_genomic_plots/"
system(paste("mkdir -p",outdir))
system("mkdir -p ../results/stackedplots")

#READING STDIN
#PosArgs <- "SARSCOV2_c_variantes_assemblies"
PosArgs <- as.character(commandArgs(trailingOnly = TRUE))
anytaxon = PosArgs[1]

#A WORD THAT THE ASSEMBLIES FILE HAS BUT THE REFERENCE DOES NOT
excluderef <- "assemblies"

#OUT_SUBDIR
system(paste("mkdir -p ",outdir,anytaxon,sep=""))

#BUILDING FEATURES AND COMPLEXITY TABLES FOR INDIVIDUAL ASSEMBLIES
# read .csv
tabla_features <- read.csv(paste("../results/GenFeatures_locations/",anytaxon,"_features_locations.csv", sep=''), header = TRUE, na.strings = "")
# Coerce Beginning and End to numeric
tabla_features$Beginning <- as.numeric(tabla_features$Beginning)
tabla_features$End <- as.numeric(tabla_features$End)
#reading complexity table
# read .csv
complexity_table <- read.csv(paste("../results/seg/",anytaxon,"_complexity.csv",sep=""),header = TRUE)
# get aa length for sequence
complexity_table$length_in_aa <- complexity_table$seg_end - complexity_table$seg_begin + 1
# get nucleotide length for sequence
complexity_table$length_in_nc <- complexity_table$length_in_aa * 3
# adjusted coordinates for seg output
complexity_table$adj_beg <- complexity_table$origin_beg + ((complexity_table$seg_begin-complexity_table$codon_start) * 3)
complexity_table$adj_end <- complexity_table$adj_beg + complexity_table$length_in_nc
# identify between high and low complexity
complexity_table$is_high <- str_detect(complexity_table$secuence, "^[:upper:]+$")
# remove duplicates
complexity_table <- complexity_table %>% distinct(viral_species, secuence, .keep_all = TRUE)
# subset to low complexity
low_complexity <- complexity_table[which(complexity_table$is_high == FALSE),] 
# plotting parameters
max_in_low_complexity <- max(low_complexity$complexity,na.rm=TRUE)
correction_axis <- 20
minimum_y_axis <- (min(low_complexity$complexity - max_in_low_complexity, na.rm=TRUE)/correction_axis)

# vector containing unique virus species
virus_species <- complexity_table$viral_species[!duplicated(complexity_table$viral_species)]
alternate_v <- as.vector(rbind(rep(-1,length(tabla_features$End)),(1))) # Only a vector to alternate positions

# EMPTY anytaxon INDEPENDENT DATAFRAMES
Anytaxon_df <- NULL
Anytaxon_df_source <- NULL
Anytaxon_df_CDSs <- NULL
Anytaxon_df_genes <- NULL
Anytaxon_df_UTRs <- NULL
Anytaxon_df_peptides <- NULL
Anytaxon_df_stemloops <- NULL
Anytaxon_df_complexity <- NULL
Anytaxon_hi_complexity <- NULL
Anytaxon_lo_complexity <- NULL
outfilename <- ""
current_genome <- ""

#FILL COUNTRY, COLLECTION AND ISOLATE
tabla_features <- tabla_features %>% fill(country)
tabla_features <- tabla_features %>% fill(collection_date)
tabla_features <- tabla_features %>% fill(isolate)

#BUILD UNIQUE VIRUS NAMES
tabla_features$Viral_species <- paste(tabla_features$Viral_species,"  ",tabla_features$country,"   ",tabla_features$collection_date," ",sep = "")
tabla_features$Viral_species <- gsub("[^a-zA-Z0-9]","_",tabla_features$Viral_species)
i = 0

#
#STACKEDPLOT
# vector to iterate in
taxids <- as.integer(complexity_table$ncbi_taxid[!duplicated(complexity_table$ncbi_taxid)])
#taxids <- as.integer(2697049)
t = taxids
# building dataframe containing all inclusive info
Anytaxon_df <- tabla_features[which(tabla_features$NCBI_taxid == t),]
Anytaxon_df <- Anytaxon_df[order(Anytaxon_df$Region_feature_class, Anytaxon_df$Beginning, Anytaxon_df$End),]
Anytaxon_df <- Anytaxon_df[!duplicated(Anytaxon_df),]
Anytaxon_df_complexity <- complexity_table[which(complexity_table$ncbi_taxid == t),]
# subsets of feature class
Anytaxon_df_source <- Anytaxon_df[which(Anytaxon_df$Region_feature_class == "source"),]
Anytaxon_df_CDSs <- Anytaxon_df[which(Anytaxon_df$Region_feature_class == "CDS"),]
Anytaxon_df_genes <- Anytaxon_df[which(Anytaxon_df$Region_feature_class == "gene"),]
Anytaxon_df_UTRs <- Anytaxon_df[which(Anytaxon_df$Region_feature_class %in% c("3_UTR","5_UTR")),]
Anytaxon_df_peptides <- Anytaxon_df[which(Anytaxon_df$Region_feature_class == "mat_peptide"),]
Anytaxon_df_stemloops <- Anytaxon_df[which(Anytaxon_df$Region_feature_class == "stem_loop"),]
# retain specific columns for each feature class
Anytaxon_df_source <- Anytaxon_df_source[c(1,2,3,4,5)]
Anytaxon_df_CDSs <- Anytaxon_df_CDSs[c(1,2,3,4,5,8)]
Anytaxon_df_genes <- Anytaxon_df_genes[c(1,2,3,4,5,8)]
Anytaxon_df_UTRs <- Anytaxon_df_UTRs[c(1,2,3,4,5)]
Anytaxon_df_peptides <- Anytaxon_df_peptides[c(1,2,3,4,5,9)]
Anytaxon_df_stemloops <- Anytaxon_df_stemloops[c(1,2,3,4,5)]
# remove duplicates
Anytaxon_df_source <- Anytaxon_df_source[!duplicated(Anytaxon_df_source), ]
Anytaxon_df_CDSs <- Anytaxon_df_CDSs[!duplicated(Anytaxon_df_CDSs), ]
Anytaxon_df_genes <- Anytaxon_df_genes[!duplicated(Anytaxon_df_genes), ]
Anytaxon_df_UTRs <- Anytaxon_df_UTRs[!duplicated(Anytaxon_df_UTRs), ]
Anytaxon_df_UTRs <- Anytaxon_df_UTRs %>% distinct(Region_feature_class, .keep_all = TRUE)
Anytaxon_df_peptides <- Anytaxon_df_peptides[!duplicated(Anytaxon_df_peptides), ]
Anytaxon_df_stemloops <- Anytaxon_df_stemloops[!duplicated(Anytaxon_df_stemloops), ]
# working with complexity
Anytaxon_df_complexity$normalized_complexity <- Anytaxon_df_complexity$complexity/max(Anytaxon_df_complexity$complexity,na.rm=TRUE)
Anytaxon_hi_complexity <- Anytaxon_df_complexity[which(Anytaxon_df_complexity$is_high == TRUE),]
Anytaxon_lo_complexity <- Anytaxon_df_complexity[which(Anytaxon_df_complexity$is_high == FALSE),]
# Fill empty lines
if(length(Anytaxon_df_stemloops$Viral_species)==0){
  Anytaxon_df_stemloops <- Anytaxon_df_stemloops %>% add_row(Viral_species = Anytaxon_df$Viral_species[1], NCBI_taxid=Anytaxon_df$NCBI_taxid[1], Region_feature_class="stem_loop", Beginning=0,End=0)
}
if(length(Anytaxon_df_source$Viral_species)==0){
  Anytaxon_df_source <- Anytaxon_df_source %>% add_row(Viral_species = Anytaxon_df$Viral_species[1], NCBI_taxid=Anytaxon_df$NCBI_taxid[1], Region_feature_class="source", Beginning=min(Anytaxon_df$Beginning,na.rm=TRUE),End=max(Anytaxon_df$End,na.rm=TRUE))
}
if(length(Anytaxon_df_CDSs$Viral_species)==0){
  Anytaxon_df_CDSs <- Anytaxon_df_CDSs %>% add_row(Viral_species = Anytaxon_df$Viral_species[1], NCBI_taxid=Anytaxon_df$NCBI_taxid[1], Region_feature_class="CDS", Beginning=0,End=0,gene="")
}
if(length(Anytaxon_df_genes$Viral_species)==0){
  Anytaxon_df_genes <- Anytaxon_df_genes %>% add_row(Viral_species = Anytaxon_df$Viral_species[1], NCBI_taxid=Anytaxon_df$NCBI_taxid[1], Region_feature_class="gene", Beginning=0,End=0,gene="")
}
if(length(Anytaxon_df_UTRs$Viral_species)==0){
  Anytaxon_df_UTRs <- Anytaxon_df_UTRs %>% add_row(Viral_species = Anytaxon_df$Viral_species[1], NCBI_taxid=Anytaxon_df$NCBI_taxid[1], Region_feature_class="", Beginning=0,End=0)
  Anytaxon_df_UTRs <- Anytaxon_df_UTRs %>% add_row(Viral_species = Anytaxon_df$Viral_species[1], NCBI_taxid=Anytaxon_df$NCBI_taxid[1], Region_feature_class="", Beginning=0,End=0)
}
if(length(Anytaxon_df_peptides$Viral_species)==0){
  Anytaxon_df_peptides <- Anytaxon_df_peptides %>% add_row(Viral_species = Anytaxon_df$Viral_species[1], NCBI_taxid=Anytaxon_df$NCBI_taxid[1], Region_feature_class="mat_peptide", Beginning=0,End=0,product="")
}
if(length(Anytaxon_df_UTRs$Viral_species)==1){
  Anytaxon_df_UTRs <- Anytaxon_df_UTRs %>% add_row(Viral_species = Anytaxon_df$Viral_species[1], NCBI_taxid=Anytaxon_df$NCBI_taxid[1], Region_feature_class="", Beginning=0,End=0)
}
if(length(Anytaxon_lo_complexity$viral_species)==0){
  Anytaxon_lo_complexity <- Anytaxon_lo_complexity %>% add_row(viral_species= Anytaxon_df$Viral_species[1], ncbi_taxid= Anytaxon_df$NCBI_taxid[1], origin_beg= 0, origin_end= 1, codon_start= 1, product= "none", seg_begin= 0, seg_end= 1, complexity= 0, window= 0, locut= 0, hicut= 0, secuence="", length_in_aa=0, length_in_nc=0,adj_beg=0,adj_end=0)
}
namesofvirus<-data.frame(unique(Anytaxon_lo_complexity$viral_species))
colnames(namesofvirus)<-c("nombres")
positions_in_y <- seq(from=0.04, to=0.15, length.out = length(namesofvirus$nombres))

sortingbyname<- sort(unique(Anytaxon_lo_complexity$viral_species))
sorting_df <- data.frame(sortingbyname,1:length(sortingbyname))
low_complexity_order <- c()

for(i in 1:length(Anytaxon_lo_complexity$viral_species)){
  Anytaxon_lo_complexity$position[i]<-which(sorting_df==Anytaxon_lo_complexity$viral_species[i])
}
for(i in 1:length(namesofvirus$nombres)){
  namesofvirus$position[i]<-which(sorting_df==namesofvirus$nombres[i])
}

f0 <- function(word,letter){
  sum(charToRaw(word)==charToRaw(letter))
}


aminoacids <- c("a","r","n","d","c","e","q","g","h","i","l","k","m","f","p","s","t","w","y","v")
aminocolor <- c("chocolate1",
                "deepskyblue4",
                "chartreuse2",
                "deeppink",
                "darkgoldenrod2",
                "deeppink4",
                "deeppink2",
                "aquamarine",
                "deepskyblue",
                "coral1",
                "coral",
                "deepskyblue2",
                "darkorange",
                "coral4",
                "coral3",
                "aquamarine4",
                "blue",
                "brown1",
                "chartreuse4",
                "brown3")
aminocoloring<-data.frame(aminoacids,aminocolor)
aminoacidcount<-c()

for(i in 1:length(Anytaxon_lo_complexity$secuence)){
  Anytaxon_lo_complexity$a[i]<-f0(Anytaxon_lo_complexity$secuence[i],"a")
  Anytaxon_lo_complexity$r[i]<-f0(Anytaxon_lo_complexity$secuence[i],"r")
  Anytaxon_lo_complexity$n[i]<-f0(Anytaxon_lo_complexity$secuence[i],"n")
  Anytaxon_lo_complexity$d[i]<-f0(Anytaxon_lo_complexity$secuence[i],"d")
  Anytaxon_lo_complexity$c[i]<-f0(Anytaxon_lo_complexity$secuence[i],"c")
  Anytaxon_lo_complexity$q[i]<-f0(Anytaxon_lo_complexity$secuence[i],"e")
  Anytaxon_lo_complexity$e[i]<-f0(Anytaxon_lo_complexity$secuence[i],"q")
  Anytaxon_lo_complexity$g[i]<-f0(Anytaxon_lo_complexity$secuence[i],"g")
  Anytaxon_lo_complexity$h[i]<-f0(Anytaxon_lo_complexity$secuence[i],"h")
  Anytaxon_lo_complexity$i[i]<-f0(Anytaxon_lo_complexity$secuence[i],"i")
  Anytaxon_lo_complexity$l[i]<-f0(Anytaxon_lo_complexity$secuence[i],"l")
  Anytaxon_lo_complexity$k[i]<-f0(Anytaxon_lo_complexity$secuence[i],"k")
  Anytaxon_lo_complexity$m[i]<-f0(Anytaxon_lo_complexity$secuence[i],"m")
  Anytaxon_lo_complexity$f[i]<-f0(Anytaxon_lo_complexity$secuence[i],"f")
  Anytaxon_lo_complexity$p[i]<-f0(Anytaxon_lo_complexity$secuence[i],"p")
  Anytaxon_lo_complexity$s[i]<-f0(Anytaxon_lo_complexity$secuence[i],"s")
  Anytaxon_lo_complexity$t[i]<-f0(Anytaxon_lo_complexity$secuence[i],"t")
  Anytaxon_lo_complexity$w[i]<-f0(Anytaxon_lo_complexity$secuence[i],"w")
  Anytaxon_lo_complexity$y[i]<-f0(Anytaxon_lo_complexity$secuence[i],"y")
  Anytaxon_lo_complexity$v[i]<-f0(Anytaxon_lo_complexity$secuence[i],"v")
}
secuencia<-""
aminocount<-c()
for(i in 1:length(Anytaxon_lo_complexity$secuence)){
  secuencia <- Anytaxon_lo_complexity$secuence[i]
  aminocount<-c(f0(secuencia,"a"),
                f0(secuencia,"r"),
                f0(secuencia,"n"),
                f0(secuencia,"d"),
                f0(secuencia,"c"),
                f0(secuencia,"e"),
                f0(secuencia,"q"),
                f0(secuencia,"g"),
                f0(secuencia,"h"),
                f0(secuencia,"i"),
                f0(secuencia,"l"),
                f0(secuencia,"k"),
                f0(secuencia,"m"),
                f0(secuencia,"f"),
                f0(secuencia,"p"),
                f0(secuencia,"s"),
                f0(secuencia,"t"),
                f0(secuencia,"w"),
                f0(secuencia,"y"),
                f0(secuencia,"v"))
  count_of_maxaa<-max(aminocount)
  maxaa<-aminoacids[which.max(aminocount)]
  #desempate
  Anytaxon_lo_complexity$maxaa[i]<-maxaa
}

for(i in 1:length(Anytaxon_lo_complexity$maxaa)){
  Anytaxon_lo_complexity$colorsito[i]<-aminocoloring$aminocolor[which(aminocoloring$aminoacids==Anytaxon_lo_complexity$maxaa[i])]
}

xlim_1=max(Anytaxon_df$End, na.rm = TRUE)
positions_in_x <- seq(from=min(Anytaxon_lo_complexity$adj_end,na.rm=TRUE), to=max(Anytaxon_lo_complexity$adj_beg,na.rm=TRUE), length.out = length(aminoacids))
mincomplex<-min(Anytaxon_lo_complexity$complexity,na.rm=TRUE)
positions_in_y

#"JUSTIFY GENOMIC PLOTS" START TO LEFT
tojustify <- list(unique(tabla_features$Viral_species))
for (i in unique(tabla_features$Viral_species)){
  terminos <- tabla_features$End[which(tabla_features$Viral_species==i)]
  tojustify[i]<-terminos
} #Recover the length of each virus

#Longest_genome <- max

stackedplot <- ggplot() + 
  scale_y_continuous(limits=c(minimum_y_axis,0.2),name="") +
  geom_rect(data=Anytaxon_lo_complexity, mapping=aes(xmin=adj_beg, xmax=adj_end, ymin=(0+(complexity-max_in_low_complexity)/correction_axis), ymax=0),fill="blue", alpha=0.07)+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))+
  geom_segment(aes(x = 0, y = 0, xend = xlim_1, yend = 0))+
  geom_segment(aes(x = 0, y = 0, xend = (1.2*xlim_1), yend = 0), color = "white",size=0.1)+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank())+
  geom_text(hjust="outward",data=namesofvirus, aes(x=max(Anytaxon_df$End), y=positions_in_y[position], label=paste(" ",gsub("[^a-zA-Z0-9]"," ",unique(nombres)))),color="black" ,size=0.5)+
  geom_text(data=Anytaxon_df, aes(x=0, y=(-0.02), label="0"),color="gray" ,size=1.5)+
  geom_text(data=Anytaxon_df, aes(x=max(End,na.rm=TRUE), y=(-0.02), label=max(End,na.rm=TRUE)),color="gray" ,size=1.5)+
  #geom_point(data=Anytaxon_lo_complexity, mapping=aes(x=adj_beg, y=positions_in_y[position]),shape=20, fill="blue", color=Anytaxon_lo_complexity$colorsito, size=0.2)+
  geom_text(data=Anytaxon_lo_complexity, aes(x=adj_beg, y=positions_in_y[position],label=maxaa ),color=Anytaxon_lo_complexity$colorsito, size=1.5)
  

#tiff("stackedplot", units="in", width=10, height=5, res=600) #tiff resolution parameters
#print(stackedplot)
#dev.off()

#### STACKEDPLOT ####
#if the taxid of the assembly is found from previous results, then rebuild reference slimcomplot
features_file_with_ref <- system(paste("grep -c ",t, "../results/GenFeatures_locations/*.csv | grep -v",excluderef,"| cut -d: -f1"), intern = TRUE)
complexity_file_with_ref <- system(paste("grep -c ",t, "../results/seg/*.csv | grep -v",excluderef,"| cut -d: -f1"), intern = TRUE)
# build two input files containing only lines with reference genome
system(paste("head -n 1",complexity_file_with_ref,"> complextmp.csv"))
system(paste("head -n 1",features_file_with_ref,"> featurestmp.csv"))
system(paste("grep",t,complexity_file_with_ref, ">> complextmp.csv"))
system(paste("grep",t,features_file_with_ref, ">> featurestmp.csv"))
# read temporary tabla_features and complexity table
tabla_featuresref <- read.csv("featurestmp.csv", header = TRUE)
tabla_featuresref$Beginning <- as.numeric(tabla_featuresref$Beginning)
tabla_featuresref$End <- as.numeric(tabla_featuresref$End)
complexity_tableref <- read.csv("complextmp.csv",header = TRUE)
system("rm *tmp.csv")
# repeat slimcomplot of reference
complexity_tableref$length_in_aa <- complexity_tableref$seg_end - complexity_tableref$seg_begin + 1
complexity_tableref$length_in_nc <- complexity_tableref$length_in_aa * 3
complexity_tableref$adj_beg <- complexity_tableref$origin_beg + ((complexity_tableref$seg_begin-complexity_tableref$codon_start) * 3)
complexity_tableref$adj_end <- complexity_tableref$adj_beg + complexity_tableref$length_in_nc
complexity_tableref$is_high <- str_detect(complexity_tableref$secuence, "^[:upper:]+$")
complexity_tableref <- complexity_tableref %>% distinct(ncbi_taxid, secuence, .keep_all = TRUE)

low_complexityref <- complexity_tableref[which(complexity_tableref$is_high == FALSE),] 
max_in_low_complexityref <- max(low_complexityref$complexity,na.rm=TRUE)
correction_axis <- 20
minimum_y_axisref <- (min(low_complexityref$complexity - max_in_low_complexityref, na.rm=TRUE)/correction_axis)
taxidsref <- complexity_tableref$ncbi_taxid[!duplicated(complexity_tableref$ncbi_taxid)] # Make a vector containing unique taxidsref
alternate_vref <- as.vector(rbind(rep(-1,length(tabla_featuresref$End)),(1))) # Only a vector to alternate positions
Anyref_dfref <- NULL
Anyref_dfref_source <- NULL
Anyref_dfref_CDSs <- NULL
Anyref_dfref_genes <- NULL
Anyref_dfref_UTRs <- NULL
Anyref_dfref_peptides <- NULL
Anyref_dfref_stemloops <- NULL
outfilename <- ""
current_genome <- ""
Anyref_dfref_complexity <- NULL
Anyref_hi_complexity <- NULL
Anyref_lo_complexity <- NULL
for(t in taxidsref){
  Anyref_dfref <- tabla_featuresref[which(tabla_featuresref$NCBI_taxid == t),]
  Anyref_dfref <- Anyref_dfref[order(Anyref_dfref$Region_feature_class, Anyref_dfref$Beginning, Anyref_dfref$End),]
  Anyref_dfref <- Anyref_dfref[!duplicated(Anyref_dfref),]
  Anyref_dfref_complexity <- complexity_tableref[which(complexity_tableref$ncbi_taxid == t),]
  # subsets of feature class
  Anyref_dfref_source <- Anyref_dfref[which(Anyref_dfref$Region_feature_class == "source"),]
  Anyref_dfref_CDSs <- Anyref_dfref[which(Anyref_dfref$Region_feature_class == "CDS"),]
  Anyref_dfref_genes <- Anyref_dfref[which(Anyref_dfref$Region_feature_class == "gene"),]
  Anyref_dfref_UTRs <- Anyref_dfref[which(Anyref_dfref$Region_feature_class %in% c("3_UTR","5_UTR")),]
  Anyref_dfref_peptides <- Anyref_dfref[which(Anyref_dfref$Region_feature_class == "mat_peptide"),]
  Anyref_dfref_stemloops <- Anyref_dfref[which(Anyref_dfref$Region_feature_class == "stem_loop"),]
  # retain specific columns for each feature class
  Anyref_dfref_source <- Anyref_dfref_source[c(1,2,3,4,5)]
  Anyref_dfref_CDSs <- Anyref_dfref_CDSs[c(1,2,3,4,5,8)]
  Anyref_dfref_genes <- Anyref_dfref_genes[c(1,2,3,4,5,8)]
  Anyref_dfref_UTRs <- Anyref_dfref_UTRs[c(1,2,3,4,5)]
  Anyref_dfref_peptides <- Anyref_dfref_peptides[c(1,2,3,4,5,9)]
  Anyref_dfref_stemloops <- Anyref_dfref_stemloops[c(1,2,3,4,5)]
  # remove duplicates
  Anyref_dfref_source <- Anyref_dfref_source[!duplicated(Anyref_dfref_source), ]
  Anyref_dfref_CDSs <- Anyref_dfref_CDSs[!duplicated(Anyref_dfref_CDSs), ]
  Anyref_dfref_genes <- Anyref_dfref_genes[!duplicated(Anyref_dfref_genes), ]
  Anyref_dfref_UTRs <- Anyref_dfref_UTRs[!duplicated(Anyref_dfref_UTRs), ]
  Anyref_dfref_UTRs <- Anyref_dfref_UTRs %>% distinct(Region_feature_class, .keep_all = TRUE)
  Anyref_dfref_peptides <- Anyref_dfref_peptides[!duplicated(Anyref_dfref_peptides), ]
  Anyref_dfref_stemloops <- Anyref_dfref_stemloops[!duplicated(Anyref_dfref_stemloops), ]
  # working with complexity
  Anyref_dfref_complexity$normalized_complexity <- Anyref_dfref_complexity$complexity/max(Anyref_dfref_complexity$complexity,na.rm=TRUE)
  Anyref_hi_complexity <- Anyref_dfref_complexity[which(Anyref_dfref_complexity$is_high == TRUE),]
  Anyref_lo_complexity <- Anyref_dfref_complexity[which(Anyref_dfref_complexity$is_high == FALSE),]
  # Fill empty lines
  if(length(Anyref_dfref_stemloops$Viral_species)==0){
    Anyref_dfref_stemloops <- Anyref_dfref_stemloops %>% add_row(Viral_species = Anyref_dfref$Viral_species[1], NCBI_taxid=Anyref_dfref$NCBI_taxid[1], Region_feature_class="stem_loop", Beginning=0,End=0)
  }
  if(length(Anyref_dfref_source$Viral_species)==0){
    Anyref_dfref_source <- Anyref_dfref_source %>% add_row(Viral_species = Anyref_dfref$Viral_species[1], NCBI_taxid=Anyref_dfref$NCBI_taxid[1], Region_feature_class="source", Beginning=min(Anyref_dfref$Beginning,na.rm=TRUE),End=max(Anyref_dfref$End,na.rm=TRUE))
  }
  if(length(Anyref_dfref_CDSs$Viral_species)==0){
    Anyref_dfref_CDSs <- Anyref_dfref_CDSs %>% add_row(Viral_species = Anyref_dfref$Viral_species[1], NCBI_taxid=Anyref_dfref$NCBI_taxid[1], Region_feature_class="CDS", Beginning=0,End=0,gene="")
  }
  if(length(Anyref_dfref_genes$Viral_species)==0){
    Anyref_dfref_genes <- Anyref_dfref_genes %>% add_row(Viral_species = Anyref_dfref$Viral_species[1], NCBI_taxid=Anyref_dfref$NCBI_taxid[1], Region_feature_class="gene", Beginning=0,End=0,gene="")
  }
  if(length(Anyref_dfref_UTRs$Viral_species)==0){
    Anyref_dfref_UTRs <- Anyref_dfref_UTRs %>% add_row(Viral_species = Anyref_dfref$Viral_species[1], NCBI_taxid=Anyref_dfref$NCBI_taxid[1], Region_feature_class="", Beginning=0,End=0)
    Anyref_dfref_UTRs <- Anyref_dfref_UTRs %>% add_row(Viral_species = Anyref_dfref$Viral_species[1], NCBI_taxid=Anyref_dfref$NCBI_taxid[1], Region_feature_class="", Beginning=0,End=0)
  }
  if(length(Anyref_dfref_peptides$Viral_species)==0){
    Anyref_dfref_peptides <- Anyref_dfref_peptides %>% add_row(Viral_species = Anyref_dfref$Viral_species[1], NCBI_taxid=Anyref_dfref$NCBI_taxid[1], Region_feature_class="mat_peptide", Beginning=0,End=0,product="")
  }
  if(length(Anyref_dfref_UTRs$Viral_species)==1){
    Anyref_dfref_UTRs <- Anyref_dfref_UTRs %>% add_row(Viral_species = Anyref_dfref$Viral_species[1], NCBI_taxid=Anyref_dfref$NCBI_taxid[1], Region_feature_class="", Beginning=0,End=0)
  }
  if(length(Anyref_lo_complexity$viral_species)==0){
    Anyref_lo_complexity <- Anyref_lo_complexity %>% add_row(viral_species= Anyref_dfref$Viral_species[1], ncbi_taxid= Anyref_dfref$NCBI_taxid[1], origin_beg= 0, origin_end= 1, codon_start= 1, product= "none", seg_begin= 0, seg_end= 1, complexity= 0, window= 0, locut= 0, hicut= 0, secuence="", length_in_aa=0, length_in_nc=0,adj_beg=0,adj_end=0)
  }
}


stackedplotref <- stackedplot + 
  geom_segment(data = Anyref_dfref_peptides , aes(x = Beginning, y = 0, xend = Beginning, yend = 0.18), color = "gray",size=0.1)+
  geom_segment(data = Anyref_dfref_peptides , aes(x = End, y = 0, xend = End, yend = 0.18),color = "gray",size=0.1)+
  geom_text(data=Anyref_dfref_peptides, aes(x=Beginning+(End-Beginning)/2, y=0.17, label=str_wrap(product, width = 8)), size=1, angle = 90)+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))+
  geom_rect(data=Anyref_dfref_genes, mapping=aes(xmin=Beginning, xmax=End, ymin=0, ymax=0.02, fill=gene), alpha=0.7)+
  geom_text(data=Anyref_dfref_genes, aes(x=Beginning+(End-Beginning)/2, y=0.01, label=str_wrap(gene, width = 10)), size=1, angle = 90)+
  geom_segment(data = Anyref_dfref_genes , aes(x = End, y = 0, xend = End, yend = 0.16),color = "gray",size=0.1)+
  theme(legend.position="none")+
  geom_segment(aes(x = 0, y = 0, xend = max(Anyref_dfref$End, na.rm = TRUE), yend = 0))+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank())+
  geom_point(data=Anyref_lo_complexity, mapping=aes(x=adj_beg, y=0.155),shape=20, fill="black", color="black", size=0.2)+
  geom_rect(aes(xmin=(min(positions_in_x))/1.5, xmax=(max(positions_in_x))*1.02, ymin=0.021, ymax=0.038), fill="gray95")+
  geom_point(aes(x=positions_in_x, y=0.03),shape=20,colour=aminocoloring$aminocolor, size=1.5)+
  geom_text(aes(x=positions_in_x, y=0.027),label=aminocoloring$aminoacids,size=1.5)

stackedplot
stackedplotref

tiff(paste("../results/stackedplots/",anytaxon,"_stackedplot.tiff",sep = ""), units="in", width=10, height=5, res=600) #tiff resolution parameters
print(stackedplotref)
dev.off()

#### heatmap ####

#building matrix
subsetcomplexity <- Anytaxon_lo_complexity[21:39]
subsetcomplexity.mx <- as.matrix(subsetcomplexity)
heatmap(subsetcomplexity.mx)


