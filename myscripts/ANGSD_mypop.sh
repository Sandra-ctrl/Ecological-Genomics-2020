myrepo="/users/s/n/snnadi/Ecological-Genomics-2020"

mkdir ${myrepo}/myresults/ANGSD

output="${myrepo}/myresults/ANGSD"

mypop="KOS"

ls /data/project_data/RS_ExomeSeq/mapping/BWA/${mypop}_*sorted.rm*.bam >${output}/${mypop}_bam.list