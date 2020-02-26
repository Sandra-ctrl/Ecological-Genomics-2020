
#!/bin/bash

cd ~/Ecological_Genomics_2020/myresults/

# i am creating a new dir to store my results

mkdir fastqcTT

for file in /data/project_data/RS_RNASeq/fastq/cleanreads/JAY_02_H*.cl.fq

do

fastqc ${file} -o fastqcTT/

done

for file in /data/project_data/RS_RNASeq/fastq/cleanreads/KAN_04_C*.cl.fq

do

fastqc ${file} -o fastqcTT/

done


