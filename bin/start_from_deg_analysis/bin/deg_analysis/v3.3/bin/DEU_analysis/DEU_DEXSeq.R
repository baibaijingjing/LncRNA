#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript


#************************section 1 usage*********************************************** 

usage=function(){
  cat("
      Descript: 
      Email:    wangmin@biomarker.com.cn
    
      Options:
            infile:  character       the path.txt file maked by bam_count.R  [forced]
            group:   character       the group information, like T01,T02_T03,T04,T05  will auto match count_list.txt [forced]
            od:      character       the output dir [forced]
            FDR:     float           the FDR cutoff to HTMLReport [default: 0.01]
        
     usage: 
          Rscript DEU_DEXSeq.R infile=path.txt group=T01,T02_T03,T04,T05 od=./
      \n"
  )
  q(status=1)
}

#************************section 2 main function*********************************************** 

DEU_DEXSeq=function(infile,group,od="./",FDR=0.01,small=F){
  library(DEXSeq)
  #--------------------------------------------
  fls=read.table(infile,as.is=T)[,1]
  gff=fls[2]
  
  #------------------split group----------------------------
  group1=strsplit(group,"_vs_")[[1]][1]
  group1=strsplit(group1,"_")[[1]]
  
  group2=strsplit(group,"_vs_")[[1]][2]
  group2=strsplit(group2,"_")[[1]]
  group=c(group1,group2)
  
  #--------------read count file and match group------------------------------
  count.fl=read.table(fls[1],as.is=T)[,1]
  
  count.flm=NULL
  for(i in 1:length(group)){
    count.flm[i]=count.fl[grepl(group[i],count.fl)]
  }
 
  #---------------------------make sample info--------------------------------
  cat("start dexseq and the sample info is :\n")
  condition=c(rep("case",length(group1)),rep("control",length(group2)))
  sampleTable = data.frame(row.names =group,condition =condition)
  sampleinfo=data.frame(file=count.flm,sampleTable)
  
  print(sampleinfo)
  cat("\n")
  
  #---------------------start dexseq--------------------------------
  dxd = DEXSeqDataSetFromHTSeq(
    count.flm,
    sampleData=sampleTable,
    flattenedfile=gff
  )
  dxd$condition <- relevel(dxd$condition, "control")
  
  #extract small data to test.  just used in test
  if(small){
    cat("extract small data to test. just used in test \n")
    dxd = dxd[geneIDs( dxd ) %in% sample(geneIDs( dxd ),100),]
  }
  
  cat("end read data, start DEU\n")
  dxr=DEXSeq(dxd,BPPARAM=MulticoreParam())
  
  cat("end DEU, start HTMLreport\n")
  DEXSeqHTML(dxr,FDR=FDR,path=paste(od,"/DEXSeqReport",sep=""),BPPARAM=MulticoreParam())
  
  d=as.data.frame(dxr)
  
  #
  write.table(d, file=paste(od,"DEU_Result_All.xls",sep=""),
              row.names = F,sep="\t",quote=F)
  
  #filter result
  d=subset(d,padj<=FDR)
  d=d[,c("groupID","featureID","log2fold_control_case","pvalue","padj")]
  colnames(d)=c("geneID","exonID","log2(FC)","pvalue","FDR")
  write.table(d, file=paste(od,"DEU_Result_Final.xls",sep=""),
              row.names = F,sep="\t",quote=F)

}

#************************section 3 call main function*********************************************** 

arg <- commandArgs(T)
#
if(length(arg)==0) usage()

# ##eval character options to R variable
arg=read.table(text=arg,sep="=",row.names=1,as.is=T)
arg=data.frame(t(arg),stringsAsFactors=F)
arg=as.list(arg)
for(i in names(arg)){arg[[i]]=type.convert(arg[[i]],as.is=T)}

#-------------------------------------------------------------
if(!all(c("infile","group","od")%in% names(arg))) {usage()}

do.call(DEU_DEXSeq,arg)
