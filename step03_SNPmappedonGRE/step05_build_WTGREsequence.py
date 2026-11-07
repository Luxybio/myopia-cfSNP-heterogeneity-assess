# Replace <LOCAL_DATA_ROOT>/<SERVER_DATA_ROOT> with your own paths.
#-*- coding: utf-8 -*-
chrnum=1
flank=30
element="intron"
allelementdir="<SERVER_DATA_ROOT>/data/20210728-data/EAS/intron/buildWTMTseq/"

saminfo=dict()
samdir="<SERVER_DATA_ROOT>/data/20210728-data/EAS/SNPonchr/flank"+str(flank)+"onchr"+str(chrnum)+".txt"
samfile=open(samdir,"r")
samline=samfile.readline()
while samline:
    samline=samline.strip("\n")
    samk=samline.split("\t")
    samkey=samk[0]+"\t"+str(samk[1])+"\t"+str(samk[3])+"\t"+str(samk[4])+"\t"+str(samk[5])
    saminfo[samkey]=1
    samline=samfile.readline()

i=1
elementdir=allelementdir+element+"onchr/"+element+"on"+str(chrnum)+".txt"
elementfile=open(elementdir,"r")
elementline=elementfile.readline()
while elementline:
    elmentline=elementline.strip("\n")
    elementinfo=elementline.split("\t")
    elementid=elementinfo[8].split(";")
    ENST=elementid[0].split(":")
    ENSE=elementid[1].split("=")
    ENSE[1]=ENSE[1].strip("\n")

    fafiledir="<SERVER_DATA_ROOT>/data/NCBIgenomeseq/Homo_GRCh38.p13_chromosome_"+str(chrnum)+".fa"
    fafile=open(fafiledir,"r")
    seqlist = fafile.read().split(">")
    for chrseq in seqlist:
        if chrseq:
            chrseq_name = chrseq.split("\n")[0]
            chrn=chrseq_name.split(" ")
            chrnc=chrn[0].split("-")
            chrseqfa = "".join(chrseq.split("\n")[1:])
            chrseqfa=chrseqfa.strip("\n")
            for rs in saminfo.keys():
                rsloc=rs.split("\t")
                if rsloc[1]=="0":
                    snploc=int(rsloc[2])-1
                    if elementinfo[6]=="+":
                        # SNP located inside an exon
                        if int(snploc)<int(elementinfo[4]) and int(snploc)>int(elementinfo[3]):
                            snponexon=int(snploc)-int(elementinfo[3])+1
                            snponexonfiledir=allelementdir+"/SNPon"+element+"loc/SNPon"+str(chrnum)+element+"loc.txt"
                            snponexonfile=open(snponexonfiledir,"a")
                            snponexonfile.write(rsloc[0]+"\t"+"0"+"\t"+ENSE[1]+"\t"+ENST[1]+"\t"+str(snponexon)+"\t"+str(snploc)+"\t"+str(rsloc[3])+"\t"+str(rsloc[2])+"\t"+str(elementinfo[3])+"\t"+str(elementinfo[4])+"\n")
                            elementseq=chrseqfa[(int(elementinfo[3])-1):int(elementinfo[4])]
                            elementnewseq=elementseq.upper()
                            snponexonseqdir=allelementdir+"/SNPon"+element+"seq/"+str(chrnum)+"_exonseq.txt"
                            snponexonseq=open(snponexonseqdir,"a")
                            snponexonseq.write(">"+ENSE[1]+"\n"+elementnewseq+"\n")
                        # SNP located at the end of an exon
                        if int(snploc)==int(elementinfo[4]):
                            snponexon=int(snploc)-int(elementinfo[3])+1
                            snponexonfiledir=allelementdir+"/SNPon"+element+"loc/SNPon"+str(chrnum)+element+"loc.txt"
                            snponexonfile=open(snponexonfiledir,"a")
                            snponexonfile.write(rsloc[0]+"\t"+"0"+"\t"+ENSE[1]+"\t"+ENST[1]+"\t"+str(snponexon)+"\t"+str(snploc)+"\t"+str(rsloc[3])+"\t"+str(rsloc[2])+"\t"+str(elementinfo[3])+"\t"+str(elementinfo[4])+"\n")
                            elementseq=chrseqfa[(int(elementinfo[3])-1):int(elementinfo[4])]
                            elementnewseq=elementseq.upper()
                            snponexonseqdir=allelementdir+"/SNPon"+element+"seq/"+str(chrnum)+"_exonseq.txt"
                            snponexonseq=open(snponexonseqdir,"a")
                            snponexonseq.write(">"+ENSE[1]+"\n"+elementnewseq+"\n")
                        if int(snploc)==int(elementinfo[3]):
                            snponexon=int(snploc)-int(elementinfo[3])+1
                            snponexonfiledir=allelementdir+"/SNPon"+element+"loc/SNPon"+str(chrnum)+element+"loc.txt"
                            snponexonfile=open(snponexonfiledir,"a")
                            snponexonfile.write(rsloc[0]+"\t"+"0"+"\t"+ENSE[1]+"\t"+ENST[1]+"\t"+str(snponexon)+"\t"+str(snploc)+"\t"+str(rsloc[3])+"\t"+str(rsloc[2])+"\t"+str(elementinfo[3])+"\t"+str(elementinfo[4])+"\n")
                            elementseq=chrseqfa[(int(elementinfo[3])-1):int(elementinfo[4])]
                            elementnewseq=elementseq.upper()
                            snponexonseqdir=allelementdir+"/SNPon"+element+"seq/"+str(chrnum)+"_exonseq.txt"
                            snponexonseq=open(snponexonseqdir,"a")
                            snponexonseq.write(">"+ENSE[1]+"\n"+elementnewseq+"\n")
                    # Exon on the negative strand: sequence handling differs
                    if elementinfo[6]=="-":
                        if int(snploc)<int(elementinfo[4]) and int(snploc)>int(elementinfo[3]):
                            snponexon=int(elementinfo[4])-int(snploc)+1
                            snponexonfiledir=allelementdir+"/SNPon"+element+"loc/SNPon"+str(chrnum)+element+"loc.txt"
                            snponexonfile=open(snponexonfiledir,"a")
                            snponexonfile.write(rsloc[0]+"\t"+"0"+"\t"+ENSE[1]+"\t"+ENST[1]+"\t"+str(snponexon)+"\t"+str(snploc)+"\t"+str(rsloc[3])+"\t"+str(rsloc[2])+"\t"+str(elementinfo[3])+"\t"+str(elementinfo[4])+"\n")
                            elementseq=chrseqfa[(int(elementinfo[3])-1):int(elementinfo[4])]
                            elementseq=elementseq.upper()
                            # Reverse complement
                            complement_AT = elementseq.replace("A","t")
                            complement_TA = complement_AT.replace("T","a")
                            complement_GC = complement_TA.replace("G","c")
                            complement_CG = complement_GC.replace("C","g")
                            complement_DNA = complement_CG.upper()
                            reverse_complement_DNA = complement_DNA[::-1]
                            snponexonseqdir=allelementdir+"/SNPon"+element+"seq/"+str(chrnum)+"_exonseq.txt"
                            snponexonseq=open(snponexonseqdir,"a")
                            snponexonseq.write(">"+ENSE[1]+"\n"+reverse_complement_DNA+"\n")
                        if int(snploc)==int(elementinfo[4]):
                            snponexon=int(elementinfo[4])-int(snploc)+1
                            snponexonfiledir=allelementdir+"/SNPon"+element+"loc/SNPon"+str(chrnum)+element+"loc.txt"
                            snponexonfile=open(snponexonfiledir,"a")
                            snponexonfile.write(rsloc[0]+"\t"+"0"+"\t"+ENSE[1]+"\t"+ENST[1]+"\t"+str(snponexon)+"\t"+str(snploc)+"\t"+str(rsloc[3])+"\t"+str(rsloc[2])+"\t"+str(elementinfo[3])+"\t"+str(elementinfo[4])+"\n")
                            elementseq=chrseqfa[(int(elementinfo[3])-1):int(elementinfo[4])]
                            elementseq=elementseq.upper()
                            # Reverse complement
                            complement_AT = elementseq.replace("A","t")
                            complement_TA = complement_AT.replace("T","a")
                            complement_GC = complement_TA.replace("G","c")
                            complement_CG = complement_GC.replace("C","g")
                            complement_DNA = complement_CG.upper()
                            reverse_complement_DNA = complement_DNA[::-1]
                            snponexonseqdir=allelementdir+"/SNPon"+element+"seq/"+str(chrnum)+"_exonseq.txt"
                            snponexonseq=open(snponexonseqdir,"a")
                            snponexonseq.write(">"+ENSE[1]+"\n"+reverse_complement_DNA+"\n")
                        if int(snploc)==int(elementinfo[3]):
                            snponexon=int(elementinfo[4])-int(snploc)+1
                            snponexonfiledir=allelementdir+"/SNPon"+element+"loc/SNPon"+str(chrnum)+element+"loc.txt"
                            snponexonfile=open(snponexonfiledir,"a")
                            snponexonfile.write(rsloc[0]+"\t"+"0"+"\t"+ENSE[1]+"\t"+ENST[1]+"\t"+str(snponexon)+"\t"+str(snploc)+"\t"+str(rsloc[3])+"\t"+str(rsloc[2])+"\t"+str(elementinfo[3])+"\t"+str(elementinfo[4])+"\n")
                            elementseq=chrseqfa[(int(elementinfo[3])-1):int(elementinfo[4])]
                            elementseq=elementseq.upper()
                            # Reverse complement
                            complement_AT = elementseq.replace("A","t")
                            complement_TA = complement_AT.replace("T","a")
                            complement_GC = complement_TA.replace("G","c")
                            complement_CG = complement_GC.replace("C","g")
                            complement_DNA = complement_CG.upper()
                            reverse_complement_DNA = complement_DNA[::-1]
                            snponexonseqdir=allelementdir+"/SNPon"+element+"seq/"+str(chrnum)+"_exonseq.txt"
                            snponexonseq=open(snponexonseqdir,"a")
                            snponexonseq.write(">"+ENSE[1]+"\n"+reverse_complement_DNA+"\n")

    elementline=elementfile.readline()

snponexonfile.close()
snponexonseq.close()
