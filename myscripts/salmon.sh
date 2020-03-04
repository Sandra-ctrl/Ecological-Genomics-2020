!#/bin/bash

cd /data/project_data/RS_RNASeq/fastq/cleanreads/

for file in JAY_02_H*.cl.fq

  do

    salmon quant -i /data/project_data/RS_RNASeq/ReferenceTranscriptome/Pabies_HC27_index \
	-l A \
	-r /data/project_data/RS_RNASeq/fastq/cleanreads/${file} \
	--validateMappings \
	-o /data/project_data/RS_RNASeq/salmon/cleanedreads/${file} 

  done 
  
  cd /data/project_data/RS_RNASeq/fastq/cleanreads/

for file in JAY_02_H*.cl.fq

  do

    salmon quant -i /data/project_data/RS_RNASeq/ReferenceTranscriptome/Pabies_HC_index \
	-l A \
	-r /data/project_data/RS_RNASeq/fastq/cleanreads/${file} \
	--validateMappings \
	-p 1 \
	--seqBias \
	-o /data/project_data/RS_RNASeq/salmon/HCmapping/${file} 

  done 
  
  
  cd /data/project_data/RS_RNASeq/fastq/cleanreads/

for file in JAY_02_H*.cl.fq

  do

    salmon quant -i /data/project_data/RS_RNASeq/ReferenceTranscriptome/Pabies_cds_index \
	-l A \
	-r /data/project_data/RS_RNASeq/fastq/cleanreads/${file} \
	--validateMappings \
	-p 1 \
	--seqBias \
	-o /data/project_data/RS_RNASeq/salmon/allmapping/${file} 

  done 