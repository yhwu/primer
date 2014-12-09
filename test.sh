## find linker at the begining of R1 reads
findprimer -f R1.fastq.gz -p GGCACATCGATTTCTGCGAGNNNNNNNNNNNNCTCCGCTTAAGGGACT -o testR1.txt 

## find primer at the beginning of R2 reads
findprimer -f R2.fastq.gz -p gaaaatctctagca -o testR2.txt

## join the table by read name
join testR1.txt testR2.txt | tr " " "\t" > R1R2.txt

## filter names of good reads; here we want R1 linkered and R2 primed
awk '$5>0 && $6<5 && $10>0 && $11<4' R1R2.txt > goodReadsR1R2.txt

## get read names and positions for R1 and R2 separately
awk '$5>0 && $6<5 && $10>0 && $11<4 {OFS="\t"; print $1,$4,$5,$6}' R1R2.txt > R1cuts.txt
awk '$5>0 && $6<5 && $10>0 && $11<4 {OFS="\t"; print $1,$9,$10,$11}' R1R2.txt > R2cuts.txt

## extract read from R1 and R2
trimprimer -f R1.fastq.gz -t R1cuts.txt -o R1.filt.fastq.gz
trimprimer -f R2.fastq.gz -t R2cuts.txt -o R2.filt.fastq.gz

## check existance of linker and primer in filtered file again
## they should give negative coordinates
findprimer -f R1.filt.fastq.gz -p GGCACATCGATTTCTGCGAGNNNNNNNNNNNNCTCCGCTTAAGGGACT -o testR1.filt.txt 
findprimer -f R2.filt.fastq.gz -p gaaaatctctagca -o testR2.filt.txt 


