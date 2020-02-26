#!/bin/bash

#this is where our output sam files are going to get converted into binary format (bam)
#then we are going to sort the bam files, remove the PCR duplicates and index them

#first, lets convert sam to bam

for f in ${output}/BWA/${mypop}*.sam

do

  out=${f/.sam/}
  sambamba-0.7.1-linux-static view -S --format=bam ${f} -o ${out}.bam
  samtools sort ${out}.bam -o ${out}.sorted.bam 
  
done

#now lets remove the PCR duplicates from our bam files

for file in ${output}/BWA/${mypop}*.sorted.bam

do

  f=${file/.sorted.bam/}
  sambamba-0.7.1-linux-static markdup -r -t 1 ${file} ${f}.sorted.rmdup.bam
  
done

#now to finish we will index our files

for file in ${output}/BWA/${mypop}*.sorted.rmdup.bam

do 

  samtools index ${file}
  
done

