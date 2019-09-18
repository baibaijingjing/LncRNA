### load library
require(edgeR)


#################################################################
### create DGE_list object
### here it will filter low count gene
#################################################################
fpkm <- function(counts=NULL, geneLength=NULL){
	# fpkm
	r <- counts
	# dim
	for ( i in 1:(dim(r)[2]) ){
		r[,i] <- 10^9 * r[,i] / sum(r[,i]) / geneLength
	}
	# return
	return(r)
}


#################################################################
### create DGE_list object
### here it will filter low count gene
#################################################################
new_DGE_list <- function(counts=NULL, condition=NULL, geneLength=NULL){
	# check NULL
	if( is.null(counts) ) stop("counts is NULL")
	if( is.null(condition) ) stop("condition is NULL")
	if( is.null(geneLength) ) stop("geneLength is NULL")

	# check length
	if( dim(counts)[2] != length(condition) )
		stop(paste("dim(counts)[2] != length(condition): ", 
			dim(counts)[2], " != ", length(condition), sep=""))
	if( dim(counts)[1] != length(geneLength) )
		stop(paste("dim(counts)[1] != length(geneLength): ", 
			dim(counts)[1], " != ", length(geneLength), sep=""))
	if( dim(counts)[2] < 2 )
		stop(paste("dim(counts)[2] < 2: ", dim(counts)[2], sep=""))

	# create DGE_list object
	de <- list(counts=as.matrix(counts), condition=condition, geneLength=geneLength,
		filterCounts=NULL, filterLength=NULL, fpkm=NULL, FDR=NULL, log2FC=NULL,
		regulated=NULL)

	# filter low expression gene
	#keep <- rowSums(cpm(de$counts)>0.2) >= min( table(condition) )
	keep <- rowSums(cpm(de$counts)>1) >= min( table(condition) )

	#keep <- rowSums(cpm(de$counts)>0.5) >= min( table(condition) )
	#keep <- rowSums(cpm(de$counts)>0) >= min( table(condition) )
	#keep <- rowSums(de$counts) >= 10  	# for debug need to remove
	#keep <- rowSums(fpkm(de$counts,de$geneLength)>1) > 0
	#keep <- rowSums(de$counts) > 0  	# for debug need to remove
	de$filterCounts <- de$counts[keep,]

	de$filterLength <- de$geneLength[keep]

	# calculate the fpkm
	de$fpkm <- fpkm(de$filterCounts, de$filterLength)
	rownames(de$fpkm) <- rownames(de$filterCounts)

	# return
	return(de)
}


#################################################################
### update DGE_list object
### here it update the dge information according to DGE analysis result
#################################################################
update_DGE_list <- function(de=NULL, FDR=NULL, log2FC=NULL){
	# check null
	# note: log2FC doesnot need to check, which is only for two condition
	if( is.null(de) ) stop("de is NULL")
	if( is.null(FDR) ) stop("FDR is NULL")

	# update
	de$FDR = FDR;
	de$log2FC = log2FC;

	# update the regulated
	if( is.null(log2FC) ){
		de$all_regulated <- rep("up", length(de$isDGE))
	}else{
		de$all_regulated <- ifelse(de$log2FC < 0, "down", "up")
	}
	# set name
	# NOTE: there may be only one DEG, so must first: ( rownames(de$fpkm) ), then [de$isDGE]
	#names(de$regulated) <- rownames(de$fpkm[de$isDGE,]) 	# NOTE: error for one DEG

	#names(de$all_regulated) <- rownames(de$fpkm)
	# update de
	de$fpkm <- de$fpkm[names(FDR),]

	de$filterCounts <- de$filterCounts[names(FDR),]
	de$filterLength <- de$filterLength[names(FDR)]
	names(de$all_regulated) <- rownames(de$fpkm)
	

	# return
	return(de)
}


#################################################################
### output DGE_list object
### here it will output the result of DGE
#################################################################
output_DGE_list <- function(de=NULL,all_file=NULL){
	# check null
	if( is.null(de) ) stop("de is NULL")
	if( is.null(all_file) ) stop("all_file is NULL")

	# for all
	if( is.null(de$log2FC) )
		df <- data.frame(id=rownames(de$fpkm), FDR=de$FDR)
	else
		df <- data.frame(id=rownames(de$fpkm),de$fpkm,FDR=de$FDR,  log2FC=de$log2FC,regulated=de$all_regulated)
		colnames(df) <- c("#ID",colnames(de$fpkm), "FDR", "log2FC","regulated")
	write.table(df, file=all_file, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}




