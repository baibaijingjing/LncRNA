#!/share/nas2/genome/biosoft/R/2.15.1/lib64/R/bin/Rscript

### include DGE plot functions
source("dge_plot_funcs.r")
source("dge_utils_funcs.r")


# usage function
usage <- function() {
	print("-------------------------------------------------------------------------------")
	print("Usage: Rscript deseq_analysis.r read_count.txt multi_factor_matrix.txt out.de")
	print("1) read_count.txt: the read count file for RNA-SEQ")
	print("2) multi_factor_matrix.txt: the multiple factor matrix")
	print("3) out.de: the output of DE")
	print("-------------------------------------------------------------------------------")
}


# abstract factor matrix
abstract_factors <- function(file) {
	# read factor
	f_mat <- read.table(file, header=TRUE)
	# check
	# return
	return(f_mat)
}


#std_two_condition_comp(f_mat, args, count_data)
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
	de <- new_DGE_list(counts=data, condition=con, geneLength=ori_data[,"geneLength"])
	print("create DEG_list object is over")
	print(paste("the number of origin gene is ", dim(de$counts)[1], sep=""))
	print(paste("the number of gene (after filtering) is ", dim(de$filterCounts)[1], sep=""))

	# for level 1
	lev1 <- NULL
	if( sum(con==levels(con)[1]) == 1 ) lev1 <- de$fpkm[,con==levels(con)[1]]
	else lev1 <- rowMeans( de$fpkm[,con==levels(con)[1]] )
	# for level 2
	lev2 <- NULL
	if( sum(con==levels(con)[2]) == 1 ) lev2 <- de$fpkm[,con==levels(con)[2]]
	else lev2 <- rowMeans( de$fpkm[,con==levels(con)[2]] )

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

	print("hello")
	print(head(res$padj))
	print(head(res))
	# update DGE_list object
	de <- update_DGE_list(de, FDR=res$padj, log2FC=res$log2FoldChange)

	# output DGE
	out <- output_DGE_list(de, all_file=paste(args[3],".final.xls",sep=""))

	# return
	return(1)
}




# get args
args <-commandArgs(TRUE)


# check args length
if( length(args) != 3 ) {
	print(args)
	usage()
	stop("the length of args != 3")
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
#ori_data <- read.delim(args[1], row.names = 1, header=TRUE)
count_data <- ori_data[ , as.vector(f_mat[,1]) ]
print("read count data is over")


# check the number of factors
if( dim(f_mat)[2] == 2 ) {
	# get levels and check
	condition <- f_mat[,2]
	print(condition)
	level <- length( levels(as.factor(condition)) )
	if( level < 2 ) {
		stop("the number of level < 2")
	} else{
		print("do std_two_condition_comp start ...")
		std_two_condition_comp(f_mat, args, count_data)
		print("do std_two_condition_comp start is over")
	}
} else if ( dim(f_mat)[2] > 2 ) {
	stop("do mutliple factor analysis is over")
} else {
	stop("factor matrix error: only one col")
}

