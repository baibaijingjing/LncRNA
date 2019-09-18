#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript

### include DGE plot functions
source("dge_plot_funcs.r")
source("dge_utils_funcs.r")


# usage function
usage <- function() {
	print("-------------------------------------------------------------------------------")
	print("Usage: Rscript ebseq_analysis.r read_count.txt S1,S2,S3 0.01 0 out.de")
	print("1) read_count.txt: the read count file for RNA-SEQ")
	print("2) S1,S2,S3: the sample id for condition")
	print("3) 0.01: FDR cutoff")
	print("4) 0: the cutoff of fold change(FC)")
	print("5) out.de: the output of DE")
	print("6) fpkm.txt: the fpkm file for RNA-SEQ")
	print("-------------------------------------------------------------------------------")
}


# only process factor without replicates
get_sample_name <- function(str=NULL) {
	# check 
	if( is.null(str) ) stop("str is NULL")
 
	# abstract factor
	s <- unlist( strsplit(str, ",") )

	# check length
	if( length(s) < 2 ) stop(paste("length(s) < 2: ", str, sep=""))

	# return
	return(s)
}


# get args
args <-commandArgs(TRUE)

# check args length
if( length(args) != 6 ) {
	print(args)
	usage()
	stop("the length of args !=6")
}



# load library
require(edgeR)
require(EBSeq)
require(ggplot2)


# abstract sample name
sample_name <- get_sample_name(args[2])


# read count data
print("read count data ...")
ori_data <- read.delim(args[1], row.names = 1, header=TRUE,check.names =F)
colnames(ori_data)<-read.delim(args[1], row.names = 1, header=F,check.names =F,stringsAsFactors=F,nrows=1)
ori_fpkm <- read.delim(args[6], row.names = 1, header=TRUE,check.names =F)
colnames(ori_fpkm)<-read.delim(args[6], row.names = 1, header=F,check.names =F,stringsAsFactors=F,nrows=1)

fpkm_data<-ori_fpkm[,sample_name]
head(ori_data)
count_data <- ori_data[,sample_name]
head(count_data)
print("read count data is over")
head(ori_fpkm)
head(fpkm_data)

# create DEG_list object
de <- new_DGE_list(counts=count_data,fpkm=fpkm_data,condition=as.factor(sample_name), geneLength=ori_data[,"geneLength"])
print("create DEG_list object is over")
print(paste("the number of origin gene is ", dim(de$counts)[1], sep=""))
print(paste("the number of gene (after filtering) is ", dim(de$filterCounts)[1], sep=""))

print(head(de$fpkm))
print (sample_name[1])
# plot corr png
#if( length(sample_name) == 2 )
#	p <- plot_corr(paste(args[5], ".cor.png", sep=""), de$fpkm[,sample_name[1]], de$fpkm[,sample_name[2]])

# calculate the FC cutoff
fc_cutoff <- ifelse( as.numeric(args[4])<=1, 0, log2(as.numeric(args[4])) )
print(fc_cutoff)


# Library size factor
Sizes=MedianNorm(de$filterCounts)

######## two level
if( length(sample_name) == 2 ) {
# gene expression estimates
# do estimate
#EBOut <- EBTest(Data=de$filterCounts, Conditions=as.factor(sample_name), sizeFactors=Sizes, maxround=5)
EBOut <- EBTest(Data=de$filterCounts, Conditions=as.factor(sample_name), sizeFactors=Sizes, maxround=5, Qtrm=0.5, QtrmCut=0)
# get pp matrix
PP <- GetPPMat(EBOut)
# get post FC
GeneFC <- PostFC(EBOut)
# print
head(GeneFC$PostFC)

# select DEG
isDGE <- (PP[,"PPDE"]>=1-as.numeric(args[3])) & (abs(log2(GeneFC$PostFC)) > fc_cutoff)
print(paste("de: ", length(isDGE), sep=""))
print(paste("PPDE = ", length(PP[,"PPDE"]), sep=""))
print(paste("(GeneFC$PostFC = ", length(GeneFC$PostFC), sep=""))

# check compare direction
C1_name <- strsplit(GeneFC$Direction, " ")[[1]][1]
log2_PostFC <- NULL
if( C1_name == sample_name[1] ){
	log2_PostFC <- -log2(GeneFC$PostFC)
} else {
	log2_PostFC <- log2(GeneFC$PostFC)
}


significant<-as.vector(GeneFC$PostFC)
a<-(PP[,"PPDE"]>=1-as.numeric(args[3])) & (log2_PostFC > fc_cutoff)
significant[a]<-"Up"
a<-(PP[,"PPDE"]>=1-as.numeric(args[3])) & (log2_PostFC < -fc_cutoff)
significant[a]<-"Down"
a<-((PP[,"PPDE"]<1-as.numeric(args[3])) | ((log2_PostFC >= -fc_cutoff) &(log2_PostFC <=fc_cutoff)))
significant[a]<-"Normal"
significant<-as.factor(significant)


# update DGE_list object
#de <- update_DGE_list(de, FDR=1-PP[,"PPDE"], log2FC=-log2(GeneFC$PostFC))
de <- update_DGE_list(de, FDR=1-PP[,"PPDE"], log2FC=log2_PostFC , isDGE=isDGE)



# plot MA and Volcano
#ma_vo <- plot_MA_Volcano(ma_file=paste(args[5], ".FC_count.png", sep=""), vo_file=paste(args[5], ".FC_FDR.png", sep=""), 
	#log10exp=log10(rowMeans(de$fpkm)), log2FC=de$log2FC, FDR=de$FDR, Significant=significant,yline=-log10(as.numeric(args[3])),xline=c(fc_cutoff,-fc_cutoff))

# output DGE
out <- output_DGE_list(de, de_file=paste(args[5],".DEG_final.xls",sep=""), 
	all_file=paste(args[5],".final.xls",sep=""), cluster_file=paste(args[5],".DEG_final.cluster",sep=""))

######## multiple level
} else {
# gene expression estimates
# do estimate
parti <- GetPatterns(as.factor(sample_name))
#parti <- parti[c(1,dim(parti)[1]),]
print(parti)
EBOut <- EBMultiTest(Data=de$filterCounts, NgVector=NULL, AllParti=parti, 
	Conditions=as.factor(sample_name), sizeFactors=Sizes, maxround=5)
# get pp matrix
PP <- GetMultiPP(EBOut)
# get post FC
#GeneFC <- PostFC(EBOut)
# print
#head(GeneFC$PostFC)

# select DEG
#isDGE <- (PP[,"PPDE"]>=1-as.numeric(args[3])) & (abs(log2(GeneFC$PostFC)) > fc_cutoff)
x_y <- matrix(c(names(PP$MAP), PP$MAP), ,2)
FDR <- 1-PP$PP[x_y]
names(FDR) <- names(PP$MAP)
pattern <- PP$MAP
isDGE <- (FDR <= as.numeric(args[3])) & (pattern != "Pattern1")

# update DGE_list object
de <- update_DGE_list(de, FDR=FDR, log2FC=NULL, isDGE=isDGE)

# output DEG
df <- data.frame(id=names(de$regulated), FDR=de$FDR[de$isDGE], Pattern=pattern[de$isDGE], regulated=de$regulated)
colnames(df) <- c("#ID", "FDR", "Pattern", "regulated")
write.table(df, file=paste(args[5], ".DEG_final.xls", sep=""), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# output DEG fpkm for cluster
df <- data.frame(Gene_lib=names(de$regulated), de$fpkm[de$isDGE,])
colnames(df) <- c("Gene_lib", colnames(de$fpkm))
write.table(df, file=paste(args[5],".DEG_final.cluster",sep=""), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# output pattern
df <- data.frame(pattern=parti)
colnames(df) <- colnames(parti)
write.table(df, file=paste(args[5],".pattern",sep=""), sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)



}




