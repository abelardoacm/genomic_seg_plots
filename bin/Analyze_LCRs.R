#LOADING LIBRARY
suppressMessages(library("tidyverse"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("viridis"))
suppressMessages(library("scales"))
suppressMessages(library("grid"))

PosArgs <- as.character(commandArgs(trailingOnly = TRUE))
#PosArgs <- "Alpha"
dataset_name = PosArgs[1]
infile <- paste("../results/seg/",dataset_name,"_complexity.csv",sep="")
seg_out <- read.table(infile, sep = ',', header = TRUE, na.strings = "", fill=TRUE,quote = "")
seg_out$is_high <- str_detect(seg_out$secuence, "^[:upper:]+$")
lowcomp <- seg_out[which(seg_out$is_high == FALSE),]



# unique genomes list
univirus <- unique(seg_out$viral_species)

# get the number of analyzed genomes
n <- length(univirus)

# get the names of proteins
proteins <- unique(lowcomp$product)

################################################################################
# subset LCR to spike
spike_LCR.db <- lowcomp[which(lowcomp$product=="surface_glycoprotein"),]

# Remove duplicated rows based on virus and sequence
spike_LCR.db <- spike_LCR.db %>% distinct(viral_species, secuence, .keep_all = TRUE)

# make vectors of positions for spike LCRs
spike1 <- c(2:11)      #spike1 | fvflvllplv
iota_exclusive <- c(252:264) #iota_exclusive | ggsssgwtagaaa
delta_kappa <- c(678:692)   #delta_kappa | srrrarsvasqsiia
delta_exclusive <- c(946:958)   #delta_exclusive | lqnvvnqnaqaln
spike4 <- c(1231:1253) #spike4 | mlccmtsccsclkgccscgscck

# label sequences in spike_LCR.db by seg_end
for (i in 1:length(spike_LCR.db$seg_end)){
  if (spike_LCR.db$seg_end[i] %in% spike1){
    spike_LCR.db$LCR_label[i] <- "spike1"
  }
  if (spike_LCR.db$seg_end[i] %in% delta_kappa){
    spike_LCR.db$LCR_label[i] <- "delta_kappa"
  }
  if (spike_LCR.db$seg_end[i] %in% delta_exclusive){
    spike_LCR.db$LCR_label[i] <- "delta_exclusive"
  }
  if (spike_LCR.db$seg_end[i] %in% spike4){
    spike_LCR.db$LCR_label[i] <- "spike4"
  }
  if (spike_LCR.db$seg_end[i] %in% iota_exclusive){
    spike_LCR.db$LCR_label[i] <- "iota_exclusive"
  }
}

# retrieve frequencies by genome of labeled LCRs
sink("spiketmp")
for (v in univirus){
  current_virus_low.db <- spike_LCR.db[which(spike_LCR.db$viral_species==v),]
  s1c <- length(which(current_virus_low.db$LCR_label=="spike1"))
  dkc <- length(which(current_virus_low.db$LCR_label=="delta_kappa"))
  dec <- length(which(current_virus_low.db$LCR_label=="delta_exclusive"))
  s4c <- length(which(current_virus_low.db$LCR_label=="spike4"))
  ioc <- length(which(current_virus_low.db$LCR_label=="iota_exclusive"))
  current_count <- paste(v,s1c,ioc,dkc,dec,s4c,sep = ",")
  print(current_count)
}
sink()

spike_LCR_count <- read.table('spiketmp', sep = ',', header = FALSE, fill=TRUE,quote = "")
system("rm spiketmp")

colnames(spike_LCR_count) <- c("virus","spike1","iota_exclusive","delta_kappa","delta_exclusive","spike4")
spike_LCR_count$virus <- gsub("^.*\"","",spike_LCR_count$virus)
spike_LCR_count$spike4 <- gsub("\"","",spike_LCR_count$spike4)

# TURN INTO NUMERIC
spike_LCR_count$spike1 <- as.numeric(spike_LCR_count$spike1)
spike_LCR_count$iota_exclusive <- as.numeric(spike_LCR_count$iota_exclusive)
spike_LCR_count$delta_kappa <- as.numeric(spike_LCR_count$delta_kappa)
spike_LCR_count$delta_exclusive <- as.numeric(spike_LCR_count$delta_exclusive)
spike_LCR_count$spike4 <- as.numeric(spike_LCR_count$spike4)

# TRANSFORM INTO 0 AND 1 EXCLUSIVELY
identifiers <- spike_LCR_count$virus
spike_LCR_count[spike_LCR_count >= 1] <- 1
spike_LCR_count$virus <- identifiers 

# BUILD SUMMARY TABLE
nvirus <- length(unique(spike_LCR_count$virus))
spike1_total <- sum(spike_LCR_count$spike1)
iota_exclusive_total <- sum(spike_LCR_count$iota_exclusive)
delta_kappa_total <- sum(spike_LCR_count$delta_kappa)
delta_exclusive_total <- sum(spike_LCR_count$delta_exclusive)
spike4_total <- sum(spike_LCR_count$spike4)

# WRITE RESULTS
#sink(paste(dataset_name,"_spike_LCR.csv",sep=""))
print("Dataset, Number of genomes, spike1 total, iota_exclusive total, delta_kappa total, delta_exclusive total, spike4 total")
print(paste(dataset_name,nvirus,spike1_total,iota_exclusive_total,delta_kappa_total,delta_exclusive_total,spike4_total,sep=","))
#sink()

#system(paste("cat ", dataset_name,"_spike_LCR.csv",sep=""))

################################################################################
# subset LCR to ORF1ab_polyprotein
ORF1ab_polyprotein_LCR.db <- lowcomp[which(lowcomp$product=="ORF1ab_polyprotein"),]

# Remove duplicated rows based on virus and sequence
ORF1ab_polyprotein_LCR.db <- ORF1ab_polyprotein_LCR.db %>% distinct(viral_species, secuence, .keep_all = TRUE)

# make vectors of positions for ORF1ab LCRs
nsp2_1 <- c(742:760) #LEGETLPTEVLTEEVVLKT
PLpro_1 <- c(926:943) #PPDEDEEEGDCEEEEFEP
PLpro_2 <- c(970:986) #QPEEEQEEDWLDDDSQQ
PLpro_3 <- c(998:1005) #QTTTIQTI
PLpro_4 <- c(1237:1255) #EEVTTTLEETKFLTENLLL
PLpro_5 <- c(2172:2183) #FFTLLLQLCTFT
PLpro_6 <- c(2546:2558) #KSKCEESSAKSAS
nsp4_1 <- c(3043:3057) #ISASIVAGGIVAIVV
nsp6_1 <- c(3583:3599) #LLLTILTSLLVLVQSTQ
nsp7_1 <- c(3869:3885) #SVVLLSVLQQLRVESSS
nsp7_2 <- c(3912:3920) #VSLLSVLLS
nsp8_1 <- c(3973:3989) #SEVVLKKLKKSLNVAKS
nsp8_2 <- c(4040:4052) #LDNDALNNIINNA


# label sequences in spike_LCR.db by seg_end
for (i in 1:length(ORF1ab_polyprotein_LCR.db$seg_end)){
  if (ORF1ab_polyprotein_LCR.db$seg_end[i] %in% nsp2_1){
    ORF1ab_polyprotein_LCR.db$LCR_label[i] <- "nsp2_1"
  }
  if (ORF1ab_polyprotein_LCR.db$seg_end[i] %in% PLpro_1){
    ORF1ab_polyprotein_LCR.db$LCR_label[i] <- "PLpro_1"
  }
  if (ORF1ab_polyprotein_LCR.db$seg_end[i] %in% PLpro_2){
    ORF1ab_polyprotein_LCR.db$LCR_label[i] <- "PLpro_2"
  }
  if (ORF1ab_polyprotein_LCR.db$seg_end[i] %in% PLpro_3){
    ORF1ab_polyprotein_LCR.db$LCR_label[i] <- "PLpro_3"
  }
  if (ORF1ab_polyprotein_LCR.db$seg_end[i] %in% PLpro_4){
    ORF1ab_polyprotein_LCR.db$LCR_label[i] <- "PLpro_4"
  }
  if (ORF1ab_polyprotein_LCR.db$seg_end[i] %in% PLpro_5){
    ORF1ab_polyprotein_LCR.db$LCR_label[i] <- "PLpro_5"
  }
  if (ORF1ab_polyprotein_LCR.db$seg_end[i] %in% PLpro_6){
    ORF1ab_polyprotein_LCR.db$LCR_label[i] <- "PLpro_6"
  }
  if (ORF1ab_polyprotein_LCR.db$seg_end[i] %in% nsp4_1){
    ORF1ab_polyprotein_LCR.db$LCR_label[i] <- "nsp4_1"
  }
  if (ORF1ab_polyprotein_LCR.db$seg_end[i] %in% nsp6_1){
    ORF1ab_polyprotein_LCR.db$LCR_label[i] <- "nsp6_1"
  }
  if (ORF1ab_polyprotein_LCR.db$seg_end[i] %in% nsp7_1){
    ORF1ab_polyprotein_LCR.db$LCR_label[i] <- "nsp7_1"
  }
  if (ORF1ab_polyprotein_LCR.db$seg_end[i] %in% nsp7_2){
    ORF1ab_polyprotein_LCR.db$LCR_label[i] <- "nsp7_2"
  }
  if (ORF1ab_polyprotein_LCR.db$seg_end[i] %in% nsp8_1){
    ORF1ab_polyprotein_LCR.db$LCR_label[i] <- "nsp8_1"
  }
  if (ORF1ab_polyprotein_LCR.db$seg_end[i] %in% nsp8_2){
    ORF1ab_polyprotein_LCR.db$LCR_label[i] <- "nsp8_2"
  }
}

# retrieve frequencies by genome of labeled LCRs
sink("ORF1ab_tmp")
for (v in univirus){
  current_virus_low.db <- ORF1ab_polyprotein_LCR.db[which(ORF1ab_polyprotein_LCR.db$viral_species==v),]
  ORF1ab_1 <- length(which(current_virus_low.db$LCR_label=="nsp2_1"))
  ORF1ab_2 <- length(which(current_virus_low.db$LCR_label=="PLpro_1"))
  ORF1ab_3 <- length(which(current_virus_low.db$LCR_label=="PLpro_2"))
  ORF1ab_4 <- length(which(current_virus_low.db$LCR_label=="PLpro_3"))
  ORF1ab_5 <- length(which(current_virus_low.db$LCR_label=="PLpro_4"))
  ORF1ab_6 <- length(which(current_virus_low.db$LCR_label=="PLpro_5"))
  ORF1ab_7 <- length(which(current_virus_low.db$LCR_label=="PLpro_6"))
  ORF1ab_8 <- length(which(current_virus_low.db$LCR_label=="nsp4_1"))
  ORF1ab_9 <- length(which(current_virus_low.db$LCR_label=="nsp6_1"))
  ORF1ab_10 <- length(which(current_virus_low.db$LCR_label=="nsp7_1"))
  ORF1ab_11 <- length(which(current_virus_low.db$LCR_label=="nsp7_2"))
  ORF1ab_12 <- length(which(current_virus_low.db$LCR_label=="nsp8_1"))
  ORF1ab_13 <- length(which(current_virus_low.db$LCR_label=="nsp8_2"))
  current_count <- paste(v,ORF1ab_1,ORF1ab_2,ORF1ab_3,ORF1ab_4,ORF1ab_5,ORF1ab_6,ORF1ab_7,ORF1ab_8,ORF1ab_9,ORF1ab_10,ORF1ab_11,ORF1ab_12,ORF1ab_13,sep = ",")
  print(current_count)
}
sink()

ORF1ab_LCR_count <- read.table('ORF1ab_tmp', sep = ',', header = FALSE, fill=TRUE,quote = "")
system("rm ORF1ab_tmp")

colnames(ORF1ab_LCR_count) <- c("virus","nsp2_1","PLpro_1","PLpro_2","PLpro_3","PLpro_4","PLpro_5","PLpro_6","nsp4_1","nsp6_1","nsp7_1","nsp7_2","nsp8_1","nsp8_2")
ORF1ab_LCR_count$virus <- gsub("^.*\"","",ORF1ab_LCR_count$virus)
ORF1ab_LCR_count$nsp8_2 <- gsub("\"","",ORF1ab_LCR_count$nsp8_2)

# TURN INTO NUMERIC
ORF1ab_LCR_count$nsp2_1 <- as.numeric(ORF1ab_LCR_count$nsp2_1)
ORF1ab_LCR_count$PLpro_1 <- as.numeric(ORF1ab_LCR_count$PLpro_1)
ORF1ab_LCR_count$PLpro_2 <- as.numeric(ORF1ab_LCR_count$PLpro_2)
ORF1ab_LCR_count$PLpro_3 <- as.numeric(ORF1ab_LCR_count$PLpro_3)
ORF1ab_LCR_count$PLpro_4 <- as.numeric(ORF1ab_LCR_count$PLpro_4)
ORF1ab_LCR_count$PLpro_5 <- as.numeric(ORF1ab_LCR_count$PLpro_5)
ORF1ab_LCR_count$PLpro_6 <- as.numeric(ORF1ab_LCR_count$PLpro_6)
ORF1ab_LCR_count$nsp4_1 <- as.numeric(ORF1ab_LCR_count$nsp4_1)
ORF1ab_LCR_count$nsp6_1 <- as.numeric(ORF1ab_LCR_count$nsp6_1)
ORF1ab_LCR_count$nsp7_1 <- as.numeric(ORF1ab_LCR_count$nsp7_1)
ORF1ab_LCR_count$nsp7_2 <- as.numeric(ORF1ab_LCR_count$nsp7_2)
ORF1ab_LCR_count$nsp8_1 <- as.numeric(ORF1ab_LCR_count$nsp8_1)
ORF1ab_LCR_count$nsp8_2 <- as.numeric(ORF1ab_LCR_count$nsp8_2)

# TRANSFORM INTO 0 AND 1 EXCLUSIVELY
identifiers <- ORF1ab_LCR_count$virus
ORF1ab_LCR_count[ORF1ab_LCR_count >= 1] <- 1
ORF1ab_LCR_count$virus <- identifiers 

# BUILD SUMMARY TABLE
nvirus <- length(unique(ORF1ab_LCR_count$virus))
nsp2_1_total <- sum(ORF1ab_LCR_count$nsp2_1)
PLpro_1_total <- sum(ORF1ab_LCR_count$PLpro_1)
PLpro_2_total <- sum(ORF1ab_LCR_count$PLpro_2)
PLpro_3_total <- sum(ORF1ab_LCR_count$PLpro_3)
PLpro_4_total <- sum(ORF1ab_LCR_count$PLpro_4)
PLpro_5_total <- sum(ORF1ab_LCR_count$PLpro_5)
PLpro_6_total <- sum(ORF1ab_LCR_count$PLpro_6)
nsp4_1_total <- sum(ORF1ab_LCR_count$nsp4_1)
nsp6_1_total <- sum(ORF1ab_LCR_count$nsp6_1)
nsp7_1_total <- sum(ORF1ab_LCR_count$nsp7_1)
nsp7_2_total <- sum(ORF1ab_LCR_count$nsp7_2)
nsp8_1_total <- sum(ORF1ab_LCR_count$nsp8_1)
nsp8_2_total <- sum(ORF1ab_LCR_count$nsp8_2)


# WRITE RESULTS
#sink(paste(dataset_name,"_ORF1ab_LCR.csv",sep=""))
print("dataset_name,nvirus,nsp2_1_total,PLpro_1_total,PLpro_2_total,PLpro_3_total,PLpro_4_total,PLpro_5_total,PLpro_6_total,nsp4_1_total,nsp6_1_total,nsp7_1_total,nsp7_2_total,nsp8_1_total,nsp8_2_total")
print(paste(dataset_name,nvirus,nsp2_1_total,PLpro_1_total,PLpro_2_total,PLpro_3_total,PLpro_4_total,PLpro_5_total,PLpro_6_total,nsp4_1_total,nsp6_1_total,nsp7_1_total,nsp7_2_total,nsp8_1_total,nsp8_2_total,sep=","))
#sink()

#system(paste("cat ", dataset_name,"_ORF1ab_LCR.csv",sep=""))

################################################################################
# subset LCR to envelope
envelope_LCR.db <- lowcomp[which(lowcomp$product=="envelope_protein"),]

# Remove duplicated rows based on virus and sequence
envelope_LCR.db <- envelope_LCR.db %>% distinct(viral_species, secuence, .keep_all = TRUE)

# make vectors of positions for envelope LCRs
envelope1 <- c(17:39)      #envelope1 | VLLFLAFVVFLLVTLAILTALRL

# label sequences in envelope_LCR.db by seg_end
for (i in 1:length(envelope_LCR.db$seg_end)){
  if (envelope_LCR.db$seg_end[i] %in% envelope1){
    envelope_LCR.db$LCR_label[i] <- "envelope1"
  }
}

# retrieve frequencies by genome of labeled LCRs
sink("envelopetmp")
for (v in univirus){
  current_virus_low.db <- envelope_LCR.db[which(envelope_LCR.db$viral_species==v),]
  env1 <- length(which(current_virus_low.db$LCR_label=="envelope1"))
  current_count <- paste(v,env1,sep = ",")
  print(current_count)
}
sink()

envelope_LCR_count <- read.table('envelopetmp', sep = ',', header = FALSE, fill=TRUE,quote = "")
system("rm envelopetmp")

colnames(envelope_LCR_count) <- c("virus","envelope1")
envelope_LCR_count$virus <- gsub("^.*\"","",envelope_LCR_count$virus)
envelope_LCR_count$envelope1 <- gsub("\"","",envelope_LCR_count$envelope1)

# TURN INTO NUMERIC
envelope_LCR_count$envelope1 <- as.numeric(envelope_LCR_count$envelope1)

# TRANSFORM INTO 0 AND 1 EXCLUSIVELY
identifiers <- envelope_LCR_count$virus
envelope_LCR_count[envelope_LCR_count >= 1] <- 1
envelope_LCR_count$virus <- identifiers 

# BUILD SUMMARY TABLE
nvirus <- length(unique(envelope_LCR_count$virus))
envelope_1_total <- sum(envelope_LCR_count$envelope1)

# WRITE RESULTS
#sink(paste(dataset_name,"_envelope_LCR.csv",sep=""))
print("Dataset, Number of genomes, envelope_1_total")
print(paste(dataset_name,nvirus,envelope_1_total,sep=","))
#sink()

#system(paste("cat ", dataset_name,"_envelope_LCR.csv",sep=""))

################################################################################
# subset LCR to ORF7a_protein
ORF7a_protein_LCR.db <- lowcomp[which(lowcomp$product=="ORF7a_protein"),]

# Remove duplicated rows based on virus and sequence
ORF7a_protein_LCR.db <- ORF7a_protein_LCR.db %>% distinct(viral_species, secuence, .keep_all = TRUE)

# make vectors of positions for envelope LCRs
ORF7a1 <- c(3:14)      #7a1 | IILFLALITLAT

# label sequences in ORF7a by seg_end
for (i in 1:length(ORF7a_protein_LCR.db$seg_end)){
  if (ORF7a_protein_LCR.db$seg_end[i] %in% ORF7a1){
    ORF7a_protein_LCR.db$LCR_label[i] <- "ORF7a_1"
  }
}

# retrieve frequencies by genome of labeled LCRs
sink("ORF7atmp")
for (v in univirus){
  current_virus_low.db <- ORF7a_protein_LCR.db[which(ORF7a_protein_LCR.db$viral_species==v),]
  ORF7a1_c <- length(which(current_virus_low.db$LCR_label=="ORF7a_1"))
  current_count <- paste(v,ORF7a1_c,sep = ",")
  print(current_count)
}
sink()

ORF7a1_LCR_count <- read.table('ORF7atmp', sep = ',', header = FALSE, fill=TRUE,quote = "")
system("rm ORF7atmp")

colnames(ORF7a1_LCR_count) <- c("virus","ORF7a_1")
ORF7a1_LCR_count$virus <- gsub("^.*\"","",ORF7a1_LCR_count$virus)
ORF7a1_LCR_count$ORF7a_1 <- gsub("\"","",ORF7a1_LCR_count$ORF7a_1)

# TURN INTO NUMERIC
ORF7a1_LCR_count$ORF7a_1 <- as.numeric(ORF7a1_LCR_count$ORF7a_1)

# TRANSFORM INTO 0 AND 1 EXCLUSIVELY
identifiers <- ORF7a1_LCR_count$virus
ORF7a1_LCR_count[ORF7a1_LCR_count >= 1] <- 1
ORF7a1_LCR_count$virus <- identifiers 

# BUILD SUMMARY TABLE
nvirus <- length(unique(ORF7a1_LCR_count$virus))
ORF7a1_total <- sum(ORF7a1_LCR_count$ORF7a_1)

# WRITE RESULTS
#sink(paste(dataset_name,"_ORF7a_1_LCR.csv",sep=""))
print("Dataset, Number of genomes, ORF7a_1_total")
print(paste(dataset_name,nvirus,ORF7a1_total,sep=","))
#sink()

#system(paste("cat ", dataset_name,"_ORF7a_1_LCR.csv",sep=""))

################################################################################
# subset LCR to ORF7b_protein
ORF7b_protein_LCR.db <- lowcomp[which(lowcomp$product %in% c("ORF7b","ORF7b_protein")),]

# Remove duplicated rows based on virus and sequence
ORF7b_protein_LCR.db <- ORF7b_protein_LCR.db %>% distinct(viral_species, secuence, .keep_all = TRUE)

# make vectors of positions for ORF7b LCRs
ORF7b1 <- c(9:30)      #7b1 | FYLCFLAFLLFLVLIMLIIFWF VOC/VOI exclusive

# label sequences in ORF7b_LCR.db by seg_end
for (i in 1:length(ORF7b_protein_LCR.db$seg_end)){
  if (ORF7b_protein_LCR.db$seg_end[i] %in% ORF7b1){
    ORF7b_protein_LCR.db$LCR_label[i] <- "ORF7b_1"
  }
}

# retrieve frequencies by genome of labeled LCRs
sink("ORF7btmp")
for (v in univirus){
  current_virus_low.db <- ORF7b_protein_LCR.db[which(ORF7b_protein_LCR.db$viral_species==v),]
  ORF7b1_c <- length(which(current_virus_low.db$LCR_label=="ORF7b_1"))
  current_count <- paste(v,ORF7b1_c,sep = ",")
  print(current_count)
}
sink()

ORF7b1_LCR_count <- read.table('ORF7btmp', sep = ',', header = FALSE, fill=TRUE,quote = "")
system("rm ORF7btmp")

colnames(ORF7b1_LCR_count) <- c("virus","ORF7b_1")
ORF7b1_LCR_count$virus <- gsub("^.*\"","",ORF7b1_LCR_count$virus)
ORF7b1_LCR_count$ORF7b_1 <- gsub("\"","",ORF7b1_LCR_count$ORF7b_1)

# TURN INTO NUMERIC
ORF7b1_LCR_count$ORF7b_1 <- as.numeric(ORF7b1_LCR_count$ORF7b_1)

# TRANSFORM INTO 0 AND 1 EXCLUSIVELY
identifiers <- ORF7b1_LCR_count$virus
ORF7b1_LCR_count[ORF7b1_LCR_count >= 1] <- 1
ORF7b1_LCR_count$virus <- identifiers 

# BUILD SUMMARY TABLE
nvirus <- length(unique(ORF7b1_LCR_count$virus))
ORF7b1_total <- sum(ORF7b1_LCR_count$ORF7b_1)

# WRITE RESULTS
#sink(paste(dataset_name,"_ORF7b_1_LCR.csv",sep=""))
print("Dataset, Number of genomes, ORF7b1_total")
print(paste(dataset_name,nvirus,ORF7b1_total,sep=","))
#sink()

#system(paste("cat ", dataset_name,"_ORF7b_1_LCR.csv",sep=""))

################################################################################
# subset LCR to nucleocapsid_phosphoprotein
nucleocapsid_LCR.db <- lowcomp[which(lowcomp$product %in% c("nucleocapsid_phosphoprotein")),]

# Remove duplicated rows based on virus and sequence
nucleocapsid_LCR.db <- nucleocapsid_LCR.db %>% distinct(viral_species, secuence, .keep_all = TRUE)

# make vectors of positions for nucleocapsid_phosphoprotein LCRs
nucleocapsid_1 <- c(175:208)#	GSRGGSQASSRSSSRSRNSSRNSTPGSSRGTSPA
nucleocapsid_2 <- c(211:230)#	AGNGGDAALALLLLDRLNQL
nucleocapsid_3 <- c(236:249)#	GKGQQQQGQTVTKK
nucleocapsid_4 <- c(361:379)#	KTFPPTEPKKDKKKKADET


# label sequences in nucleocapsid_phosphoprotein by seg_end
for (i in 1:length(nucleocapsid_LCR.db$seg_end)){
  if (nucleocapsid_LCR.db$seg_end[i] %in% nucleocapsid_1){
    nucleocapsid_LCR.db$LCR_label[i] <- "nucleocapsid_1"
  }
  if (nucleocapsid_LCR.db$seg_end[i] %in% nucleocapsid_2){
    nucleocapsid_LCR.db$LCR_label[i] <- "nucleocapsid_2"
  }
  if (nucleocapsid_LCR.db$seg_end[i] %in% nucleocapsid_3){
    nucleocapsid_LCR.db$LCR_label[i] <- "nucleocapsid_3"
  }
  if (nucleocapsid_LCR.db$seg_end[i] %in% nucleocapsid_4){
    nucleocapsid_LCR.db$LCR_label[i] <- "nucleocapsid_4"
  }
}

# retrieve frequencies by genome of labeled LCRs
sink("nucleocapsidtmp")
for (v in univirus){
  current_virus_low.db <- nucleocapsid_LCR.db[which(nucleocapsid_LCR.db$viral_species==v),]
  nucleo1_c <- length(which(current_virus_low.db$LCR_label=="nucleocapsid_1"))
  nucleo2_c <- length(which(current_virus_low.db$LCR_label=="nucleocapsid_2"))
  nucleo3_c <- length(which(current_virus_low.db$LCR_label=="nucleocapsid_3"))
  nucleo4_c <- length(which(current_virus_low.db$LCR_label=="nucleocapsid_4"))
  current_count <- paste(v,nucleo1_c,nucleo2_c,nucleo3_c,nucleo4_c,sep = ",")
  print(current_count)
}
sink()

nucleocapsid_LCR_count <- read.table('nucleocapsidtmp', sep = ',', header = FALSE, fill=TRUE,quote = "")
system("rm nucleocapsidtmp")

colnames(nucleocapsid_LCR_count) <- c("virus","nucleocapsid_1","nucleocapsid_2","nucleocapsid_3","nucleocapsid_4")
nucleocapsid_LCR_count$virus <- gsub("^.*\"","",nucleocapsid_LCR_count$virus)
nucleocapsid_LCR_count$nucleocapsid_4 <- gsub("\"","",nucleocapsid_LCR_count$nucleocapsid_4)

# TURN INTO NUMERIC
nucleocapsid_LCR_count$nucleocapsid_1 <- as.numeric(nucleocapsid_LCR_count$nucleocapsid_1)
nucleocapsid_LCR_count$nucleocapsid_2 <- as.numeric(nucleocapsid_LCR_count$nucleocapsid_2)
nucleocapsid_LCR_count$nucleocapsid_3 <- as.numeric(nucleocapsid_LCR_count$nucleocapsid_3)
nucleocapsid_LCR_count$nucleocapsid_4 <- as.numeric(nucleocapsid_LCR_count$nucleocapsid_4)

# TRANSFORM INTO 0 AND 1 EXCLUSIVELY
identifiers <- nucleocapsid_LCR_count$virus
nucleocapsid_LCR_count[nucleocapsid_LCR_count >= 1] <- 1
nucleocapsid_LCR_count$virus <- identifiers 

# BUILD SUMMARY TABLE
nvirus <- length(unique(nucleocapsid_LCR_count$virus))

nucleocapsid_1_total <- sum(nucleocapsid_LCR_count$nucleocapsid_1)
nucleocapsid_2_total <- sum(nucleocapsid_LCR_count$nucleocapsid_2)
nucleocapsid_3_total <- sum(nucleocapsid_LCR_count$nucleocapsid_3)
nucleocapsid_4_total <- sum(nucleocapsid_LCR_count$nucleocapsid_4)

# WRITE RESULTS
#sink(paste(dataset_name,"_nucleocapsid_LCR.csv",sep=""))
print("Dataset, Number of genomes, nucleocapsid_1_total,nucleocapsid_2_total,nucleocapsid_3_total,nucleocapsid_4_total")
print(paste(dataset_name,nvirus,nucleocapsid_1_total,nucleocapsid_2_total,nucleocapsid_3_total,nucleocapsid_4_total,sep=","))
#sink()

#system(paste("cat ", dataset_name,"_nucleocapsid_LCR.csv",sep=""))

# BUILD RESULTS
Resultados <- as.data.frame(t(data.frame(nsp2_1_total,
                         PLpro_1_total,
                         PLpro_2_total,
                         PLpro_3_total,
                         PLpro_4_total,
                         PLpro_5_total,
                         PLpro_6_total,
                         nsp4_1_total,
                         nsp6_1_total,
                         nsp7_1_total,
                         nsp7_2_total,
                         nsp8_1_total,
                         nsp8_2_total,
                         spike1_total,
                         iota_exclusive_total,
                         delta_kappa_total,
                         delta_exclusive_total,
                         spike4_total,
                         envelope_1_total,
                         ORF7a1_total,
                         ORF7b1_total,
                         nucleocapsid_1_total,
                         nucleocapsid_2_total,
                         nucleocapsid_3_total,
                         nucleocapsid_4_total)))
colnames(Resultados)[1]<-"Total_count"
rownames(Resultados) <- gsub("_total","",rownames(Resultados))
rownames(Resultados) <- gsub("_"," ",rownames(Resultados))
Resultados$Percentage <- (Resultados$Total_count / nvirus)*100
Resultados <- tibble::rownames_to_column(Resultados, "LCR")

#Turn 'treatment' column into a character vector
Resultados$LCR <- as.character(Resultados$LCR)
Resultados$LCR <- factor(Resultados$LCR, levels=unique(Resultados$LCR))

for (i in 1:length(Resultados$Percentage)){
  if (Resultados$Percentage[i] == 0){
    Resultados$Frequency[i] <- "non-existent"
  }
  if (Resultados$Percentage[i] > 0){
    Resultados$Frequency[i] <- "extremely rare"
  }
  if (Resultados$Percentage[i] > 0.5){
    Resultados$Frequency[i] <- "very rare"
  }
  if (Resultados$Percentage[i] > 1){
    Resultados$Frequency[i] <- "rare"
  }
  if (Resultados$Percentage[i] > 10){
    Resultados$Frequency[i] <- "infrequent"
  }
  if (Resultados$Percentage[i] > 30){
    Resultados$Frequency[i] <- "moderately prevalent"
  }
  if (Resultados$Percentage[i] > 60){
    Resultados$Frequency[i] <- "prevalent"
  }
  if (Resultados$Percentage[i] > 70){
    Resultados$Frequency[i] <- "very prevalent"
  }
  if (Resultados$Percentage[i] > 90){
    Resultados$Frequency[i] <- "highly prevalent"
  }
  if (Resultados$Percentage[i] > 99){
    Resultados$Frequency[i] <- "extremely prevalent"
  }
  if (Resultados$Percentage[i] == 100){
    Resultados$Frequency[i] <- "allways present"
  }
}

Resultados$Frequency <- as.character(Resultados$Frequency)
Resultados$Frequency <- factor(Resultados$Frequency, levels=rev(c("non-existent", "extremely rare", "very rare","rare","infrequent","moderately prevalent","prevalent","very prevalent","highly prevalent","extremely prevalent","allways present")))
enetext = textGrob(paste("n = ", nvirus, sep = ""))

Resultados$Percentage <- Resultados$Percentage %>% round(., 1)

p <- ggplot(Resultados, aes(x=LCR, y=Percentage, fill=Frequency))+
  geom_bar(stat="identity")+
  scale_fill_viridis_d() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(0, 100)+
  geom_text(aes(label=paste(Percentage,"%",sep="")), vjust=-0.3, size=3.5)+
  labs(y=paste('Percentage (n = ',nvirus,')',sep=""))

system("mkdir -p ../results/LCR_counts")
tiff(paste("../results/LCR_counts/",dataset_name,"_LCRcount.tiff",sep=""), units="in", width=18, height=6, res=300) #tiff resolution parameters
p
dev.off()

# LCR PRESENCE MATRIX
spike_LCR_count$virus <- NULL 
envelope_LCR_count$virus <- NULL 
ORF7a1_LCR_count$virus <- NULL 
ORF7b1_LCR_count$virus <- NULL 
nucleocapsid_LCR_count$virus <- NULL
LCR_presencia <- cbind(ORF1ab_LCR_count,spike_LCR_count,envelope_LCR_count,ORF7a1_LCR_count,ORF7b1_LCR_count,nucleocapsid_LCR_count)

write_csv(LCR_presencia,paste("../results/LCR_counts/",dataset_name,"_LCR_matrix.csv",sep=""))

system(paste("cat ../results/LCR_counts/",dataset_name,"_LCR_matrix.csv",sep=""))
