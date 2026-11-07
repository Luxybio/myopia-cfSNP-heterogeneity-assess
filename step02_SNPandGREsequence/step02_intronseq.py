# Replace <LOCAL_DATA_ROOT>/<SERVER_DATA_ROOT> with your own paths.
import re
file=open("<LOCAL_DATA_ROOT>/DNA-element-data/test.txt","r",encoding="utf-8")
# Read FASTA format and concatenate multi-line sequences into a single line
loc=[]
seq_list=file.read().split(">")
for seq in seq_list :
    if seq :
        seq_name=seq.split("\n")[0]
        seq_fa="".join(seq.split("\n")[1:])
        id=seq_name.split("|")
        # Exon start/end
        loc_start=id[8].split(";")
        loc_start.sort()
        loc_end=id[9].split(";")
        loc_end.sort()
        loc_str=loc_start+loc_end
        # Get the position of the first exon on the chromosome; used when converting to sequence length
        loc_1=min(loc_str)
        j=0
        for i in range(0 , len(loc_start)-1):
            j=j+1
            # Intron start/end
            intron_start=int(loc_end[i]) + 1
            intron_end=int(loc_start[j]) - 1
            # Extract the sequence
            all_seq=str(seq_fa)
            jie_start_index=int(intron_start) - int(loc_1)
            jie_end_index=int(intron_end) + 1 - int(loc_1)
            jie = all_seq[jie_start_index:jie_end_index]
            print(id[0])
            print (str(intron_start)+"-"+str(intron_end))
            print (str(jie_start_index)+"-"+str(jie_end_index))
            print(jie)
            out=open("<LOCAL_DATA_ROOT>/DNA-element-data/intron_test.txt","a",encoding="utf-8")
            out.write(id[0]+"_"+id[1]+"_"+id[2]+"_"+str(intron_start)+"-"+str(intron_end)+str(jie_start_index)+"-"+str(jie_end_index)+"\n")
            out.write(jie+"\n")
