# Replace <LOCAL_DATA_ROOT>/<SERVER_DATA_ROOT> with your own paths.
import re
file=open("<LOCAL_DATA_ROOT>/seek-data/ENSMBL-Biomart/Homo_sapiens.GRCh38.102.chromosome.Y.gff3.txt","r",encoding="utf-8")
exon={}
intron={}
for line in file:
    if line[0] == "#":
        continue
    else:
        line=line.strip()
        str=line.split("\t")
        # Match lines containing exons
        if re.search("exon",str[2],re.IGNORECASE):
            info=str[8].split(";")
            ENST=info[0].split(":")
            ENSE=info[5].split("=")
            id=str[0]+"|"+ENST[1]+"|"+ENSE[1]
            e_loc=str[3]+"-"+str[4]
            # Use chr|ENST|ENSE as key, exon start-end as value
            exon[id]=e_loc
            intron_id=ENST[1]
            intron_value=ENSE[1]+":"+str[3]+"-"+str[4]
            if intron_id in intron:
                intron[intron_id]=intron[intron_id]+"|"+intron_value
            else:
                intron[intron_id]=intron_value

for keys in intron.keys():
    value=intron[keys]
    ENSE_loc=value.split("|")
    start=[]
    end=[]
    for element in ENSE_loc:
        num=element.split(":")
        num_loc=num[1].split("-")
        start.append(num_loc[0])
        end.append(num_loc[1])
    start=";".join(start)
    end=";".join(end)
    intron_gap=keys+"|"+start+"|"+end
    print (intron_gap)
