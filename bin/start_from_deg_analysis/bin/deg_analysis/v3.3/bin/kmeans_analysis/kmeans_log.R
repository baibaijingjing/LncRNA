library(pheatmap)
library('getopt');

spec = matrix(c(
  'help'    , 'h', 0, "logical",
  'infile'  , 'i', 1, "character",
  'outdir'  , 'o', 1, "character",
  'cluster' , 'c', 1, "integer",
  'show_xlab' ,'a', 0, "logical",
  'rownum'  , 'r', 1, "integer",
  'colnum'  , 'n', 1, "interger",
  'alpha'  , 'alp', 1, "double"
), byrow=TRUE, ncol=4);
opt = getopt(spec);


# define usage function
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("
      Usage example: 
      Rscript  kmeans_linux.R --infile test.txt  --outdir ./
      Rscript  kmeans_linux.R --infile test.txt  --outdir ./  --cluster 10
      
      Options: 
      --help  	 -h 	NULL 		    get this help
      --infile 		character 	    the input file [forced]
      --outdir 		character 	    the dirname for output graph and file [forced]
      --rownum    character       the number of cluster for each row [auto selection]
      --colnum    character       the number of column needed [auto selection]
      --cluster 	integer 	      the number of cluster(K).[default: 6]
      --show_xlab	logical 	     show xlab [optional, default: T]
      --alpha     float          the alpha of colour  [default: 0.1]
      \n")
  q(status=1);
}
##############################################################################



##############################################################################################
x_y=function(xy){
  #---------------------------------------------------------------
  #computer the best rownum and colnum  for the number of k
  #---------------------------------------------------------------
  x=floor(xy^0.5)
  if(x==xy^0.5){
    y=x
  }else{
    y=x+1
  }
  repeat{
    if(x*y>=xy){break}
    x=x+1
  }
  c(x,y)
}

col.alpha <- function(col, alpha=1){ 
  #--------------------------------------------------------
  #reduce the alpha value (level of transparency) of colours 
  #---------------------------------------------------------
  if(missing(col)) stop("Please provide a vector of colours.") 
  apply(sapply(col, col2rgb)/255, 2,  
        function(x)rgb(x[1], x[2], x[3], alpha=alpha)
  )   
} 

#*******************************************************************************************
#end function


#----------------------------------------------start command--------------------------------
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$infile) )	{ print_usage(spec) }
if ( is.null(opt$outdir) )	{ print_usage(spec) }
if ( is.null(opt$norm) )  { opt$norm="auto" }
if ( is.null(opt$cluster) )	{ opt$cluster=6}
if ( is.null(opt$rownum) )	{opt$rownum=x_y(opt$cluster)[1]}
if ( is.null(opt$colnum) )	{opt$colnum=x_y(opt$cluster)[2]}
if ( is.null(opt$show_xlab ) )	{ opt$show_xlab =T }
if ( is.null(opt$alpha ) )  { opt$alpha =0.2}

data<-read.delim(opt$infile,sep="\t",row.names=1)  ## file 
data2=data
data=log2(data+1)

kmean=kmeans(data,opt$cluster,iter.max=1000,nstar=10)  ####cluster number 

cluster=kmean$cluster
write.table(kmean$centers,file=paste(opt$outdir,"kmeans_centers.txt",sep="/"),sep="\t",quote=F,row=F)
res_cluster=data.frame(id=rownames(data2),cluster=cluster,data2)
write.table(res_cluster,file=paste(opt$outdir,"kmeans_cluster.txt",sep="/"),sep="\t",quote=F,row=F)

#png(file=paste(opt$outdir,"kmeans_centers_heatmap.png",sep="/"),width =105, height = 76.8,unit="mm",res=600)
#pheatmap(kmean$centers,cluster_col=F,border_color=NA,color=topo.colors(200),show_colnames=F,fontsize=15)
#dev.off()

cluster<-as.matrix(cluster)
png(file=paste(opt$outdir,"k-means.png",sep="/"), width = opt$colnum*33, height = opt$rownum*30,unit="mm",res=600)
par(mfrow=c(opt$rownum,opt$colnum),mar=c(1.2, 2, 1.2, 0.4) + 0.1,mgp=c(0.6,0.05,0),tck=0.02, cex.axis=0.5,cex.lab=0.5,
    cex.main=0.7,font.lab=1,font.axis=1,font.main=1,col.axis="gray20")
cut<-col.alpha(rainbow(opt$cluster),opt$alpha)
for (i in 1:opt$cluster) {
  x<-cluster[,1]==i
  part<-data[x,]
  len<-dim(part)[1]
  part<-rbind(part,kmean$centers[i,])
  part<-rbind(part,kmean$centers[i,])
  part<-rbind(part,kmean$centers[i,])
  matplot(t(part),type='l',lty=1,col=c(rep(cut[i],len),rep("black",3)),
          axes=F,main=paste("Cluster",i,sep=""),xlab="",ylab=expression(log[2](FPKM+1)),ylim=range(data))
  
  box(lwd=0.5,col="gray20")
  axis(2,lwd=0.5,col="gray20",las=1)
  if (opt$show_xlab) {axis(1,1:dim(data)[2],labels=colnames(data),lwd=0.5,col="gray20",mgp=c(0,-0.2,0))}
}
dev.off()


#for (i in 1:opt$cluster) {
#  x<-cluster[,1]==i
#  part<-data[x,]
#  #part2<-data2[x,]
#  #write.table(part2,file=paste(opt$outdir,"/Cluster",i,'.txt',sep=""),sep="\t",quote=F)
#  png(file=paste(opt$outdir,"/Cluster",i,'_heatmap.png',sep=""),width = 1024, height = 768,)
#  pheatmap(part,cluster_col=F,border_color=NA,fontsize=15,color=topo.colors(200),show_colnames=F,show_rownames=F,main=paste("Cluster",i,sep=""))
#  dev.off()
#}