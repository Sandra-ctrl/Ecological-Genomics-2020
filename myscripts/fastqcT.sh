
#!/bin/bash

cd ~/Ecological_Genomics_2020/myresults/

# i am creating a new dir to store my results

mkdir fastqcT

#for file in /data/project_data/RS_RNASeq/fastq/JAY_02_H*fastq.gz

#do

#fastqc ${file} -o fastqcT/

#done

for file in /data/project_data/RS_RNASeq/fastq/KAN_04_C*fastq.gz

do

fastqc ${file} -o fastqcT/

done


