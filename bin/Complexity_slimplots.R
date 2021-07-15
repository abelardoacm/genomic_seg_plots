

library("tidyverse")
library("ggrepel")

outdir <- "../results/Complexity_genomic_plots/"
system(paste("mkdir -p",outdir))
# reading input
#PosArgs <- "Nidovirales"
PosArgs <- as.character(commandArgs(trailingOnly = TRUE))
anytaxon = PosArgs[1]
#Out_subdir
system(paste("mkdir -p ",outdir,anytaxon,sep=""))

tabla_features <- read.csv(paste("../results/GenFeatures_locations/",anytaxon,"_features_locations.csv", sep=''), header = TRUE)
tabla_features$Beginning <- as.numeric(tabla_features$Beginning)
tabla_features$End <- as.numeric(tabla_features$End)

#reading complexity table
complexity_table <- read.csv(paste("../results/seg/",anytaxon,"_complexity.csv",sep=""),header = TRUE)
complexity_table$length_in_aa <- complexity_table$seg_end - complexity_table$seg_begin + 1
complexity_table$length_in_nc <- complexity_table$length_in_aa * 3
complexity_table$adj_beg <- complexity_table$origin_beg + ((complexity_table$seg_begin-complexity_table$codon_start) * 3)
complexity_table$adj_end <- complexity_table$adj_beg + complexity_table$length_in_nc
complexity_table$is_high <- str_detect(complexity_table$secuence, "^[:upper:]+$")
complexity_table <- complexity_table %>% distinct(ncbi_taxid, secuence, .keep_all = TRUE)
low_complexity <- complexity_table[which(complexity_table$is_high == FALSE),] 
max_in_low_complexity <- max(low_complexity$complexity,na.rm=TRUE)
correction_axis <- 20
minimum_y_axis <- (min(low_complexity$complexity - max_in_low_complexity, na.rm=TRUE)/correction_axis)


#############################################################
taxids <- complexity_table$ncbi_taxid[!duplicated(complexity_table$ncbi_taxid)] # Make a vector containing unique taxids
alternate_v <- as.vector(rbind(rep(-1,length(tabla_features$End)),(1))) # Only a vector to alternate positions
Anytaxon_df <- NULL
# subsets of feature class
Anytaxon_df_source <- NULL
Anytaxon_df_CDSs <- NULL
Anytaxon_df_genes <- NULL
Anytaxon_df_UTRs <- NULL
Anytaxon_df_peptides <- NULL
Anytaxon_df_stemloops <- NULL
outfilename <- ""
current_genome <- ""
Anytaxon_df_complexity <- NULL
Anytaxon_hi_complexity <- NULL
Anytaxon_lo_complexity <- NULL


for(t in taxids){
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
  outfilename = paste(Anytaxon_df$Viral_species[1],"id",Anytaxon_df$NCBI_taxid[1],"slimplot")
  outfilename <- str_replace_all(outfilename, "[^[:alnum:]]", "_")
  outfilename <- paste(outfilename,".tiff",sep="")
  tiff(paste("../results/Complexity_genomic_plots/",anytaxon,"/",outfilename,sep=""), units="in", width=10, height=1.8, res=600) #tiff resolution parameters
  current_genome <- ggplot() + 
    scale_x_continuous(name=paste(Anytaxon_df$Viral_species[1])) + 
    scale_y_continuous(limits=c(minimum_y_axis,0.15),name="") +
    geom_segment(data = Anytaxon_df_peptides , aes(x = Beginning, y = 0, xend = Beginning, yend = 0.15), color = "gray",size=0.1)+
    geom_segment(data = Anytaxon_df_peptides , aes(x = End, y = 0, xend = End, yend = 0.15),color = "gray",size=0.1)+
    geom_rect(data=Anytaxon_lo_complexity, mapping=aes(xmin=adj_beg, xmax=adj_end, ymin=(0+(complexity-max_in_low_complexity)/correction_axis), ymax=0),fill="blue", alpha=0.3)+
    geom_text(data=Anytaxon_df_peptides, aes(x=Beginning+(End-Beginning)/2, y=0.1, label=str_wrap(product, width = 8)), size=1, angle = 90)+
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))+
    geom_rect(data=Anytaxon_df_genes, mapping=aes(xmin=Beginning, xmax=End, ymin=0, ymax=0.07, fill=gene), alpha=0.7)+
    geom_text(data=Anytaxon_df_genes, aes(x=Beginning+(End-Beginning)/2, y=0.035, label=str_wrap(gene, width = 10)), size=1, angle = 90)+
    theme(legend.position="none")+
    geom_segment(aes(x = 0, y = 0, xend = max(Anytaxon_df$End, na.rm = TRUE), yend = 0))+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank())+
    geom_text(vjust="inward",hjust="inward",data=Anytaxon_df, aes(x=0, y=minimum_y_axis, label=Anytaxon_df$Viral_species[1]),color="black" ,size=2)+
    geom_text(data=Anytaxon_df, aes(x=0, y=(-0.02), label="0"),color="gray" ,size=1.5)+
    geom_text(data=Anytaxon_df, aes(x=max(End,na.rm=TRUE), y=(-0.02), label=max(End,na.rm=TRUE)),color="gray" ,size=1.5)+
    geom_point(data=Anytaxon_lo_complexity, mapping=aes(x=adj_beg, y=0.07),shape=20, fill="blue", color="darkred", size=0.5)
  print(current_genome)
  dev.off()
}


