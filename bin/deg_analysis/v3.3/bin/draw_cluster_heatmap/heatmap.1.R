#!/share/nas2/genome/biosoft/R/2.15.1/lib64/R/bin/Rscript
#####################################################################
# Copyright 2015, BMK
#
# Author:tengh <tengh@biomarker.com.cn>
#
# Function: draw genomewide cytosine coverage distribution map
#
# Modify date: 20150819
# Note: delete group label
#	reset opt$color="#263C8B,#4E74A6,#BDBF78,#BFA524"
#####################################################################

# load library
library('getopt');
#opt<-data.frame(infile="E:/R_workplace/20150626heatmap/T1_T2_vs_T3_T4.DEG.final.cluster",groupfile="E:/R_workplace/20150626heatmap/groupfile.heatmap")
.sourcePath<-"/share/nas1/tengh/research/Rsource/"

#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'help' , 'h', 0, "logical",
  'infile' , 'i', 1, "character",
  'groupfile' , 'G', 2, "character",
  'outfile','o',2,"character",
  'cell.width' , 'w', 2, "double",
  'cell.height','e',2,"double",
  'title','t',2,"character",
  'width','W',2,"integer",
  'height','H',2,"integer",
  'size','s',2,"double",
  'rowname','R',2,"logical",
  'colname','C',2,"logical",
  'color','c',2,"character" ,
  'zero','z',2,"double",
  'log','l',2,"character",
  'scale','S',2,"character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

#遗传图与基因组的共线性分析
# define usage function
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("
      Usage example: 
      1) Rscript heatmap.R --infile in.heatmap --outfile heatmap --color BrBG
      2) Rscript heatmap.R --infile in.heatmap --outfile heatmap --groupfile group.heatmap --title heatmap --size 10 --rownames F
      3) Rscript heatmap.R --infile in.heatmap --outfile heatmap --title heatmap --size 10 --cell.width 7 --cell.height 7



      
      Options: 
      --help          -h  NULL        get this help
      --infile        -i  character   the tab delimited input file saving numeric matrix of the values to be plotted.[forced]
      --outfile       -o  character   file path where to save the picture. Filetype is decided by the extension in the path. [optional,heatmap in current working directory]
      --groupfile     -G  character   the tab delimited input file saving data frame that specifies the annotations shown on top side of the heatmap [optional, default:NA]
      --cell.width    -w  double      individual cell width in points[optional, default: 7]
      --cell.height   -e  double      individual cell height in points[optional, default: 7]
      --size          -s  double      base fontsize for the plot[optional, default: 10]
      --width         -W  double      manual option for determining the output file width in pixel.[optional, default: NA]
      --heigth        -H  double      manual option for determining the output file height in pixel.[optional, default:NA]
      --title         -t  character   a title for the plot[optional, default: ]
      --rowname       -R  logical     boolean specifying if row names are be shown.[optional, default: TRUE]
      --colname       -C  logical     boolean specifying if column names are be shown.[optional, default:NA]
      --color         -c  character   choose the colour set(redgreen BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral)or set colour splited by , .[optional, default: BrBG]
      --zero          -z  double      Set the minima vlaue: set mat values less than minima to minima.[optional, default:1]
      --log           -l  character   a logarithmic log scale is in use.[optional, default:log2]
      --scale         -S  character   character indicating if the values should be centered and scaled in either the row direction or the column direction, or none..[optional, default:none]
      \n")
  q(status=1);
}

#if(file.exists(paste(.sourcePath,"heatmap/pheatmap2.r",sep="")))source(paste(.sourcePath,"heatmap/pheatmap2.r",sep=""))else stop(paste(.sourcePath,"heatmap/pheatmap2.r does not exist!",sep=""))

# if help was asked for print a friendly message
# and exit with a non-zero error code
if (!is.null(opt$help)) { print_usage(spec) }


# check non-null args
if ( is.null(opt$infile) )  { print_usage(spec) }else {opt$infile<-gsub("\\\\",replacement = "/",x = opt$infile)}

#set some reasonable defaults for the options that are needed,
#but were not specified.
if ( is.null(opt$groupfile) )  { opt$groupfile=NA }else {opt$groupfile<-gsub("\\\\",replacement = "/",x = opt$groupfile)}

if( is.null(opt$outfile))opt$outfile="heatmap"
if(is.null(opt$title))opt$title=""

if(is.null(opt$width)){
  opt$width=NA  
}else if(!(is.numeric(opt$width)&&opt$width>0)){
  stop("Parameter Error：outfile width must be positive integer")  
}else{
  opt$width=opt$width/500
}

if(is.null(opt$height)){
  opt$height=NA
}else if(!(is.numeric(opt$height)&&opt$height>0)){
  stop("Parameter Error：outfile height must be positive integer")  
}else{
  opt$height=opt$height/500
}


if(is.null(opt$cell.width)){
  opt$cell.width=ifelse(is.na(opt$width),7,NA)
}else if(!(is.numeric(opt$cell.width)&&opt$cell.width>0)){
  stop("Parameter Error：cell width must be positive integer")
}

if(is.null(opt$cell.height)){
  opt$cell.height=ifelse(is.na(opt$height),7,NA )  
}else if(!(is.numeric(opt$cell.height)&&opt$cell.height>0)){
  stop("Parameter Error：cell height must be positive integer")
}

if(is.null(opt$rowname))opt$rowname=T
if(is.null(opt$colname))opt$colname=T
#if(is.null(opt$color))opt$color="RdYlGn"
if(is.null(opt$color))opt$color="#263C8B,#4E74A6,#BDBF78,#BFA524"
if(is.null(opt$zero))opt$zero=1
if(is.null(opt$size))opt$size=10
if(is.null(opt$log))opt$log="log2"
if(is.null(opt$scale))opt$scale="none"

##import data
rawdat<-read.table(opt$infile,head=T,sep="\t",comment.char = "")
#rawdat<-read.table(as.vector(opt$infile),header=T,sep="\t",comment.char = "")
rownames(rawdat)<-as.matrix(rawdat)[,1]
rawdat<-as.matrix(rawdat[,-1])
rawdat<-rawdat+opt$zero

if(opt$log=="log2"){  
  rawdat<-log2(rawdat)
}else if(opt$log=="log10"){
  rawdat<-log10(rawdat)
}else if(is.na(opt$log)){
  rawdat=rawdat
}else{
  stop("Paramter error: a logarithmic scale parameter log can only be NA log10 or log2!")
}


#
if(is.na(opt$groupfile)){
  anColor = NA
  colGroup =NA
  heat.dat=rawdat
}else{
  groupdat<-read.table(as.vector(opt$groupfile),header=F,sep="\t",comment.char = "")
  group<-as.vector(groupdat[,2])
  names(group)<-as.vector(groupdat[,1])
  if(sum(!is.element(names(group),colnames(rawdat)))>0){
    stop(paste(c("the following samples in group file not exist:",setdiff(names(group),colnames(rawdat)),"please check your groupfile!"),sep="\n"))
  }
  if(sum(!is.element(colnames(rawdat),names(group)))>0){
    warning(paste(c("the following samples in infile will not be ploted:",setdiff(names(group),colnames(rawdat))),sep="\n"))
  }
  #多类样品热图添加分类条
  heat.dat<-rawdat[,names(group)]
  colGroup<-data.frame(Group=group)
  colGroup$Group= factor(colGroup$Group, levels = c(unique(group), "other"))
  row.names(colGroup)<-names(group)#设置样品颜色类
  gColor<-c( "#7FC97F","#BEAED4","#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666", "#B3E2CD","#FDCDAC","#CBD5E8","#F4CAE4","#E6F5C9","#FFF2AE","#F1E2CC","#CCCCCC")
  gColor=gColor[1:length(unique(group))]
  names(gColor)<-unique(group)
  anColor<-list(Group=gColor)
}


if(length(opt$color)==1&&is.element(opt$color,c("BrBG","PiYG","PRGn","PuOr","RdBu","RdGy","RdYlBu","RdYlGn","Spectral"))){  
  require(RColorBrewer)
  hmColors=colorRampPalette(rev(brewer.pal(n = 7, name = opt$color)))(100)
}else if(length(opt$color)==1&&(opt$color=="redgreen")){
  library(gplots) 
  message(paste("color=",opt$color,sep=""))
  hmColors=redgreen(255)
}else{
  hmColors<-strsplit(opt$color,split = ",")[[1]]
  hmColors=colorRampPalette(hmColors)(256)    
}
hl<-hclust(dist(heat.dat))
capture.output(str(as.dendrogram(hclust(dist(heat.dat)))),file =paste(c(opt$outfile,".txt"),collapse =""))
#message(c("width",opt$width,"height",opt$height))
pheatmap2(filename =paste(c(opt$outfile,".png"),collapse =""),width = opt$width,height = opt$height,mat=heat.dat,cellwidth=opt$cell.width,color = hmColors,cellheight=opt$cell.height,main=opt$title,cluster_rows=T,cluster_cols=T,annotation_col = colGroup,annotation = colGroup,annotation_colors = anColor,border_color=NA,fontsize=opt$size,col=hmColors,show_rownames=opt$rowname,show_colnames=opt$colname,fontsize_col=ifelse(is.na(opt$cell.width),opt$size,min(opt$size,opt$cell.width)),fontsize_row=ifelse(is.na(opt$cell.height),opt$size,min(opt$size,opt$cell.height)),scale=opt$scale)
dev.off()
pheatmap2(filename =paste(c(opt$outfile,".pdf"),collapse =""),width = opt$width,height = opt$height,mat=heat.dat,cellwidth=opt$cell.width,color = hmColors,cellheight=opt$cell.height,main=opt$title,cluster_rows=T,cluster_cols=T,annotation_col = colGroup,annotation = colGroup,annotation_colors = anColor,border_color=NA,fontsize=opt$size,col=hmColors,show_rownames=opt$rowname,show_colnames=opt$colname,fontsize_col=ifelse(is.na(opt$cell.width),opt$size,min(opt$size,opt$cell.width)),fontsize_row=ifelse(is.na(opt$cell.height),opt$size,min(opt$size,opt$cell.height)),scale=opt$scale)
dev.off()
