#!/usr/bin/env python3
# Replace <LOCAL_DATA_ROOT>/<SERVER_DATA_ROOT> with your own paths.
# Extract 30 bp upstream and downstream of the SNP position for TF binding site analysis

file=open("<SERVER_DATA_ROOT>/data/new-DNA-SNP-function-11/SNPonDNA-alle-loc.txt","r")

line=file.readline()
while line:
   line=line.strip("\n")
   info=line.split("\t")
   ensr=info[0]
   loc=int(info[4])
   # Start/end positions
   startseq=int(loc)-1-30
   endseq=int(loc)-1+30+1

   wsfiledir="<SERVER_DATA_ROOT>/data/new-DNA-SNP-function-11/WT-regulatoryelement-seq/"+ensr+".txt"
   msfiledir="<SERVER_DATA_ROOT>/data/new-DNA-SNP-function-11/MT-regulatoryelement-seq/"+ensr+".txt"
   # Extract 30 bp upstream/downstream of the wild-type
   wsfile=open(wsfiledir,"r")
   wseqlist = wsfile.read().split(">")
   for wsseq in wseqlist:
    if wsseq:
        wsseq_name = wsseq.split("\n")[0]
        wsseqfa = "".join(wsseq.split("\n")[1:])

        wnewseq=wsseqfa[int(startseq):int(endseq)]
        woutdir="<SERVER_DATA_ROOT>/data/new-DNA-SNP-function-11/regulatoryelement-localseq/WT-regulatoryelement-seq/"+ensr+".txt"
        wout=open(woutdir,"w")
        wout.write(">"+ensr+"\n"+wnewseq+"\n")

   # Extract 30 bp upstream/downstream of the mutant
   msfile=open(msfiledir,"r")
   mseqlist = msfile.read().split(">")
   for msseq in mseqlist:
    if msseq:
        msseq_name = msseq.split("\n")[0]
        msseqfa = "".join(msseq.split("\n")[1:])

        mnewseq=msseqfa[int(startseq):int(endseq)]
        moutdir="<SERVER_DATA_ROOT>/data/new-DNA-SNP-function-11/regulatoryelement-localseq/MT-regulatoryelement-seq/"+ensr+".txt"
        mout=open(moutdir,"w")
        mout.write(">"+ensr+"\n"+mnewseq+"\n")

   line=file.readline()
