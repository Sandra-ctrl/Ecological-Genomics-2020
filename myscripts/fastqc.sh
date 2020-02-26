
#!/bin/bash

cd ~/Ecological_Genomics_2020/myresults/

# i am creating a new dir to store my results

mkdir fastqc

for file in /data/project_data/RS_ExomeSeq/fastq/edge_fastq/KOS*fastq.gz

do

fastqc ${file} -o fastqc/

done


