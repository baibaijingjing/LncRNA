#!/share/nas2/genome/biosoft/R/2.15.1/lib64/R/bin/Rscript

### include DGE plot functions
source("dge_plot_funcs.r")
source("dge_utils_funcs.r")


# usage function
usage <- function() {
	print("-------------------------------------------------------------------------------")
	print("Usage: Rscript ebseq_analysis.r read_count.txt S1,S2,S3 out.de")
	print("1) read_count.txt: the read count file for RNA-SEQ")
	print("2) S1,S2,S3: the sample id for condition")
	print("3) out.de: the output of DE")
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
if( length(args) != 3 ) {
	print(args)
	usage()
	stop("the length of args != 3")
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
#ori_data <- read.delim(args[1], row.names = 1, header=TRUE)
head(ori_data)
count_data <- ori_data[,sample_name]
head(count_data)
print("read count data is over")


# create DEG_list object
de <- new_DGE_list(counts=count_data, condition=as.factor(sample_name), geneLength=ori_data[,"geneLength"])

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

# check compare direction
C1_name <- strsplit(GeneFC$Direction, " ")[[1]][1]
log2_PostFC <- NULL
if( C1_name == sample_name[1] ){
	log2_PostFC <- -log2(GeneFC$PostFC)
} else {
	log2_PostFC <- log2(GeneFC$PostFC)
}

# update DGE_list object
#de <- update_DGE_list(de, FDR=1-PP[,"PPDE"], log2FC=-log2(GeneFC$PostFC))
de <- update_DGE_list(de, FDR=1-PP[,"PPDE"], log2FC=log2_PostFC)



# output DGE
out <- output_DGE_list(de, all_file=paste(args[3],".final.tmp.xls",sep=""))

} 



