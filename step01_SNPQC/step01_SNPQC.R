# Replace <LOCAL_DATA_ROOT>/<SERVER_DATA_ROOT> with your own paths.
rm(list=ls())
strGsub = function(x,REF,ALT){
  x = gsub(pattern = "0", replacement = x[REF], x = x)
  x = gsub(pattern = "1", replacement = x[ALT], x = x)
  x = gsub(pattern = "\\|", replacement = "/", x = x)
  x = gsub(pattern = "\\.", replacement = "", x = x)
  x
}
library(genetics)
calcHWP = function(x,start,end){
  g = genotype(x[start:end],allow.partial.missing=TRUE)
  snphweRes = HWE.chisq(g,simulate.p.value=FALSE)
  snphwePvalue = snphweRes$p.value
  snphwePvalue
}
callratio = function(x,start,end){
  g = genotype(x[start:end],allow.partial.missing=TRUE)
  ratio = (length(g)-sum(is.na(g)))/length(g)
  ratio
}
genotyperatio = function(x){
  g = genotype(x,allow.partial.missing=TRUE)
  ratio = (length(g)-sum(is.na(g)))/length(g)
  ratio
}

path = "<LOCAL_DATA_ROOT>/20201231_data_analysis/EUR/"
outputpath = "<LOCAL_DATA_ROOT>/20201231_data_analysis/test_out/"
filenames = list.files(path=path, pattern="ALL.chr1_GRCh38.genotypes.20170504.vcf")
i=1
files = paste(path,filenames[i],sep="")
print(files)
outputfiles = paste(outputpath,filenames[i],sep="")
data = read.table(files,sep="\t",header=TRUE,stringsAsFactors = FALSE,comment.char = "",colClasses="character")
snpinfo = data[,c(4:5,10:dim(data)[2])]
snpinfo = apply(snpinfo,1,strGsub,REF="REF",ALT="ALT")
snpinfo = t(snpinfo)
# call ratio > 75%
call_ratio = apply(snpinfo,1,callratio,start=3,end = dim(snpinfo)[2])
snpinfo = snpinfo[which(call_ratio>0.75),]
data = data[which(call_ratio>0.75),]
# genotype > 75%
genotype_ratio = apply(snpinfo[,3:dim(snpinfo)[2]],2,genotyperatio)
snpinfo = snpinfo[,c(1:2,which(genotype_ratio>0.75)+2)]
data = data[,c(1:9,which(genotype_ratio>0.75)+9)]
# H-W P-value > 0.01
index=c()
for (j in (1:length(snpinfo[,1])))
{
  if(max(table(snpinfo[j,]))!=(length(snpinfo[1,])-2))
  {
    index=c(index,j)
  }
}
snpinfo=snpinfo[index,]
data=data[index,]
HWPvalue = apply(snpinfo,1,calcHWP,start=3,end = dim(snpinfo)[2])
snpinfo = snpinfo[which(HWPvalue>0.01),]
data = data[which(HWPvalue>0.01),]
print (filenames[i])
write.table(data,file=outputfiles,sep="\t",col.names=TRUE,row.names = FALSE,quote=FALSE)
