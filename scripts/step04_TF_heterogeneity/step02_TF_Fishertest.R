# Replace <LOCAL_DATA_ROOT>/<SERVER_DATA_ROOT> with your own paths.

rm(list=ls())
# Fisher's exact test
 dir="<LOCAL_DATA_ROOT>/27-result-buchong/TF/TF-xishu-WT-MT-Fisher/"
file=list.files(dir)
for (i in 1:length(file)){
  filedir=paste0(dir,file[i])
  data=read.table(filedir,sep="\t")

  data <- data.frame(in.sum=c(data[,1],data[,2]),no.sum=c(data[,3],data[,4]))

  rownames(data)=c("WT","MT")

  res=fisher.test(data,alternative="two.sided")
  resp=res$p.value
  resci=res$conf.int
  print(file[i])
  print(resp)

  strend=paste(file[i],resp,collapse="\t")
  print(strend)
  # Add ratio of WT and MT bound TF counts
  rfile=paste0("<LOCAL_DATA_ROOT>/27-result-buchong/TF/WT-MT-TF-ratio/",file[i])
  rdata=read.table(rfile)
  rnum=rdata[1,2]

  strend1=paste(strend,rnum,collapse="\t")

  write.table(strend1,"<LOCAL_DATA_ROOT>/27-result-buchong/TF/WT-MT-fisher-pvalue.txt",append =T,
              quote=F,col.names = F,row.names=F,sep="\t")
}


# Visualization
rm(list=ls())
library(ggpubr)
library(ggbreak)
library(ggrepel)
library(ggplot2)

pfile=read.table("<LOCAL_DATA_ROOT>/27-result-buchong/TF/WT-MT-fisher-pvalue-phe.txt",sep="\t")
pfile[,5]=-log10(pfile[,2]) # -log10(p-value)

colnames(pfile)=c("id","p","ratio","phe","log")

p=ggscatter(pfile, x = "ratio", y = "log",palette = "nature",shape="phe",mean.point = FALSE,star.plot = FALSE,color="red")
p4=p+geom_hline(aes(yintercept=1.30))
p5=p4+geom_text_repel(label=pfile$id)
p1=p+scale_y_break(c(5,10))
p2=p1+geom_text_repel(label=pfile$id)
p3=p2+geom_hline(aes(yintercept=1.30))
