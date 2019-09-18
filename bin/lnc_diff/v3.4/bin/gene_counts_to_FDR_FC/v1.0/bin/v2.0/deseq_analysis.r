#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript

### include DGE plot functions
source("dge_plot_funcs.r")
source("dge_utils_funcs.r")


# usage function
usage <- function() {
	print("-------------------------------------------------------------------------------")
	print("Usage: Rscript deseq_analysis.r read_count.txt multi_factor_matrix.txt 0.01 0 out.de")
	print("1) read_count.txt: the read count file for RNA-SEQ")
	print("2) multi_factor_matrix.txt: the multiple factor matrix")
	print("3) 0.01: FDR cutoff")
	print("4) 0: the cutoff of log 2 fold change(log2FC)")
	print("5) out.de: the output of DE")
	print("6) fpkm.txt: the fpkm file for RNA-SEQ")
	print("-------------------------------------------------------------------------------")
}


# abstract factor matrix
abstract_factors <- function(file) {
	# read factor
	f_mat <- read.table(file, header=TRUE,check.names =F)
	# check
	# return
	return(f_mat)
}



# compare two condition
# One factor with two levels
# note: must one have biological replicates
std_two_condition_comp <- function(f_mat, args, data) {
	# get levels and check
	condition <- f_mat[,2]
	print(condition)
	level <- length( levels(as.factor(condition)) )
	if( level != 2 ) {
		stop("the level of factor for std_two_condition_comp must == 2")
	}

	# create DEG_list object
	con <- as.factor(condition)
	print(con)
	de <- new_DGE_list(counts=data,fpkm=fpkm_data, condition=con, geneLength=ori_data[,"geneLength"])
	print("create DEG_list object is over")
	print(paste("the number of origin gene is ", dim(de$counts)[1], sep=""))
	print(paste("the number of gene (after filtering) is ", dim(de$filterCounts)[1], sep=""))

	# plot corr png
	# for level 1
	lev1 <- NULL
	if( sum(con==levels(con)[1]) == 1 ) lev1 <- de$fpkm[,con==levels(con)[1]]
	else lev1 <- rowMeans( de$fpkm[,con==levels(con)[1]] )
	# for level 2
	lev2 <- NULL
	if( sum(con==levels(con)[2]) == 1 ) lev2 <- de$fpkm[,con==levels(con)[2]]
	else lev2 <- rowMeans( de$fpkm[,con==levels(con)[2]] )
	# plot
	#p <- plot_corr(paste(args[5], ".cor.png", sep=""), lev1, lev2)

	# create count data set
	cds <- newCountDataSet(countData = de$filterCounts, conditions = condition)
	dim(cds)

	## Normalization
	cds <- estimateSizeFactors(cds)
	sizeFactors(cds)
	## Estimate Dispersion
	#cds <- estimateDispersions(cds)
	cds <- estimateDispersions(cds, fitType="local") 	# can for miRNA or small gene number
	
	## Differential Expression
	# Identify those features that are direntially expressed in the two groups.
	print("Differential Expression ... ")
	c1 <- names( table(condition) )[1]
	c2 <- names( table(condition) )[2]
	res <- nbinomTest(cds, c1, c2)
	print("Differential Expression is over ")

	# calculate the FC cutoff
	fc_cutoff <- ifelse( as.numeric(args[4])<=1, 0, log2(as.numeric(args[4])) )
	print(fc_cutoff)

	# select DEG
	isDGE <- (res$padj<as.numeric(args[3])) & (abs(res$log2FoldChange)>fc_cutoff) 

	significant<-as.vector(res$padj)
	a<-res$padj<as.numeric(args[3]) & res$log2FoldChange > fc_cutoff
	significant[(res$padj<as.numeric(args[3])) & (res$log2FoldChange > fc_cutoff)]<-"Up"
	significant[(res$padj<as.numeric(args[3])) & (res$log2FoldChange < -fc_cutoff)]<-"Down"
	significant[(res$padj>=as.numeric(args[3]))|((res$log2FoldChange <= fc_cutoff) & (res$log2FoldChange >= -fc_cutoff)) ]<-"Normal"
	significant<-as.factor(significant)
	print("hello")
	print(head(res$padj))
	print(head(res))
	# update DGE_list object
	de <- update_DGE_list(de,isDGE=isDGE,FDR=res$padj,log2FC=res$log2FoldChange)

	# plot MA and Volcano
	#ma_vo <- plot_MA_Volcano(ma_file=paste(args[5], ".FC_count.png", sep=""), vo_file=paste(args[5], ".FC_FDR.png", sep=""),log10exp=log10(rowMeans(de$fpkm)), log2FC=de$log2FC, FDR=de$FDR, Significant=significant,yline=-log10(as.numeric(args[3])),xline=c(fc_cutoff,-fc_cutoff))
	

	# output DGE
	out <- output_DGE_list(de, de_file=paste(args[5],".DEG_tmp.xls",sep=""),
		all_file=paste(args[5],".tmp.xls",sep=""), cluster_file=paste(args[5],".DEG_final.cluster",sep=""))

	# return
	return(length(isDGE))
}


# compare multiple condition
# One factor with multiple levels
# note: must one have biological replicates
std_multiple_level_comp <- function(f_mat, args, data) {
	# get levels and check
	condition <- f_mat[,2]
	print(condition)
	level <- length( levels(as.factor(condition)) )
	if( level <= 2 ) {
		stop("the level of factor for std_multiple_level_comp must > 2")
	}

	# create DEG_list object
	con <- as.factor(condition)
	print(con)
	de <- new_DGE_list(counts=data,fpkm=fpkm_data, condition=con, geneLength=ori_data[,"geneLength"])
	print("create DEG_list object is over")
	print(paste("the number of origin gene is ", dim(de$counts)[1], sep=""))
	print(paste("the number of gene (after filtering) is ", dim(de$filterCounts)[1], sep=""))

	# create count data set
	cds <- newCountDataSet(countData=de$filterCounts, conditions = condition)
	dim(cds)

	## Normalization
	cds <- estimateSizeFactors(cds)
	sizeFactors(cds)
	## Estimate Dispersion
	cds <- estimateDispersions(cds, method="pooled") 	# diff
	
	## Differential Expression
	# Identify those features that are direntially expressed in the two groups.
	print("Differential Expression ... ")
	fit1 <- fitNbinomGLMs( cds, count ~ condition )
	fit0 <- fitNbinomGLMs( cds, count ~ 1 )
	pvalsGLM <- nbinomGLMTest( fit1, fit0 )
	padjGLM = p.adjust( pvalsGLM, method="BH" )
	print( head( rownames(padjGLM) ) )
	print( head( padjGLM ) )
	print("Differential Expression is over ")
	
	# calculate the FC cutoff
	fc_cutoff <- ifelse( as.numeric(args[4])<=1, 0, log2(as.numeric(args[4])) )
	print(fc_cutoff)

	# select DEG
	isDGE <- (padjGLM<as.numeric(args[3]))


	# update DGE_list object
	de <- update_DGE_list(de, FDR=padjGLM, log2FC=NULL, isDGE=isDGE)

	# plot MA and Volcano
	#ma_vo <- plot_MA_Volcano(ma_file=paste(args[5], ".FC_count.png", sep=""), vo_file=paste(args[5], ".FC_FDR.png", sep=""),
	#	log10exp=log10(rowMeans(de$fpkm)), log2FC=de$log2FC, FDR=de$FDR, Significant=isDGE)

	# output DGE
	out <- output_DGE_list(de, de_file=paste(args[5],".DEG_tmp.xls",sep=""),
		all_file=paste(args[5],".tmp.xls",sep=""), cluster_file=paste(args[5],".DEG_final.cluster",sep=""))

	# return
	return(length(isDGE))
}


# get args
args <-commandArgs(TRUE)


# check args length
if( length(args) != 6 ) {
	print(args)
	usage()
	stop("the length of args != 6")
}


# abstract multiple factors matrix
print("abstract_factors is start")
f_mat <- abstract_factors(args[2])
print("abstract_factors is over")


# load library
require(edgeR)
require(DESeq)
require(ggplot2)


# read Count data
print("read count data ...")
ori_data <- read.delim(args[1], row.names = 1, header=TRUE,check.names =F)
colnames(ori_data)<-read.delim(args[1], row.names = 1, header=F,check.names =F,stringsAsFactors=F,nrows=1)
count_data <- ori_data[ , as.vector(f_mat[,1]) ]

print("read count data is over")

print("fpkm data ...")
ori_fpkm <- read.delim(args[6], row.names = 1, header=TRUE,check.names =F)
colnames(ori_fpkm)<-read.delim(args[6], row.names = 1, header=F,check.names =F,stringsAsFactors=F,nrows=1)

fpkm_data <- ori_fpkm[ , as.vector(f_mat[,1]) ]
print("fpkm data is over")

# check the number of factors
if( dim(f_mat)[2] == 2 ) {
	# get levels and check
	condition <- f_mat[,2]
	print(condition)
	level <- length( levels(as.factor(condition)) )
	if( level < 2 ) {
		stop("the number of level < 2")
	} else if ( level == 2 ) {
		print("do std_two_condition_comp start ...")
		std_two_condition_comp(f_mat, args, count_data)
		print("do std_two_condition_comp start is over")
	}else {
		print("do std_multiple_level_comp start ...")
		std_multiple_level_comp(f_mat, args, count_data)
		print("do std_multiple_level_comp start is over")
	}
} else if ( dim(f_mat)[2] > 2 ) {
	stop("do mutliple factor analysis is over")
} else {
	stop("factor matrix error: only one col")
}

