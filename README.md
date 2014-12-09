findprimer
=====

This  program checks if primer or linker exists in a read straight from Illumina sequencing machine, assuming that the linker or primer is at the beginning of a read. The program works as follows. Given a primer such as gaaaatctctagca, or a linker such as GGCACATCGATTTCTGCGAGNNNNNNNNNNNNCTCCGCTTAAGGGACT, the program calculates the edit distance between the primer(linker) and the beginning part of a read, and output the coordinates where the reads should be cut. Here are examples of output: 
```
@M03249:8:000000000-ABY6R:1:1101:19266:2177     1:N:0:0 GGCACATCGATTTCTGCGAGNNNNNNNNNNNNCTCCGCTTAAGGGACT        0
       47      0
@M03249:8:000000000-ABY6R:1:1101:11955:2181     1:N:0:0 GGCACATCGATTTCTGCGAGNNNNNNNNNNNNCTCCGCTTAAGGGACT        0
       48      1
@M03249:8:000000000-ABY6R:1:1101:16499:2188     1:N:0:0 GGCACATCGATTTCTGCGAGNNNNNNNNNNNNCTCCGCTTAAGGGACT        0
       47      0
@M03249:8:000000000-ABY6R:1:1101:13346:2194     1:N:0:0 GGCACATCGATTTCTGCGAGNNNNNNNNNNNNCTCCGCTTAAGGGACT        0
       47      0
@M03249:8:000000000-ABY6R:1:1101:20494:2279     1:N:0:0 GGCACATCGATTTCTGCGAGNNNNNNNNNNNNCTCCGCTTAAGGGACT        0
       46      1
```
where 46,47,48's are 0-based last match points, and the last numbers are edit distances. One should filter the reads based on the edit distances. Normally, primer/linker length divided by 10 is a good choice.

trimprimer
=====

Given sequence names and cut points, this program selects the reads and cuts them. Note that the sequence names should be in the same order as in the fastq file, as the program only does sequential search once looking for the reads. Examples of input are listed below. 

```
@M03249:8:000000000-ABY6R:1:1101:15413:1732	0	47	0
@M03249:8:000000000-ABY6R:1:1101:16481:1751	0	48	0
@M03249:8:000000000-ABY6R:1:1101:17183:1757	0	47	0
@M03249:8:000000000-ABY6R:1:1101:18023:1770	0	46	0
```

#### A typical work flow for Illumina miseq run

- 1. demultiplex with idemp.
```
idemp 
```

- 2. Find inker at the begining of R1 reads, primer in R2 reads
```
findprimer -f R1.fastq.gz -p GGCACATCGATTTCTGCGAGNNNNNNNNNNNNCTCCGCTTAAGGGACT -o testR1.txt 
findprimer -f R2.fastq.gz -p gaaaatctctagca -o testR2.txt
```

- 3. Join the table by read name
```
join testR1.txt testR2.txt | tr " " "\t" > R1R2.txt
```

- 4. Get read names and cut positions for R1 and R2 separately
```
awk '$5>0 && $6<5 && $10>0 && $11<4 {OFS="\t"; print $1,$4,$5,$6}' R1R2.txt > R1cuts.txt
awk '$5>0 && $6<5 && $10>0 && $11<4 {OFS="\t"; print $1,$9,$10,$11}' R1R2.txt > R2cuts.txt
```

- 5. Extract read from R1 and R2
```
trimprimer -f R1.fastq.gz -t R1cuts.txt -o R1.filt.fastq.gz
trimprimer -f R2.fastq.gz -t R2cuts.txt -o R2.filt.fastq.gz
```

- 6. Check existance of linker and primer in filtered file again; they should give negative coordinates.
```
findprimer -f R1.filt.fastq.gz -p GGCACATCGATTTCTGCGAGNNNNNNNNNNNNCTCCGCTTAAGGGACT -o testR1.filt.txt 
findprimer -f R2.filt.fastq.gz -p gaaaatctctagca -o testR2.filt.txt 
```

- 7. Map reads.


#### Installation
```
git clone https://github.com/yhwu/primer
cd primer
make
make test
```


#### Usage
```
[yhwu@local primer]$ ./findprimer
Usage:
   findprimer -f fastq -p primer -m n -o outFile

Options:
   fastq    fastq file
   primer   primer or linker sequence
   n        allowed base mismatches, optional, default=1+primer/20
   outFile  output folder, optional, default=.

Output: rows of the following columns
   sequence_name
   sequence_name_comment
   primer/linker
   primer/linker_length
   primer/linker_on_sequence
   edit distance
   primer_start	#0 based
   primer_end	#0 based

Note:
   1. primer is converted to upper cases.
   2. sequence is not converted.
   3. N, ?, . in primer matches all.
   4. N in sequence does not matche any.

[yhwu@local primer]$ ./trimprimer
Usage:
   trimprimer -f fastq -t trimfile -o outFile

Options:
   fastq    fastq file
   trimfile primer or linker sequence
   outFile  output folder, optional, default=.

Note: the trimfile should contain rows with the following field
   sequence_name  start with @
   start          ignored, always cut from 0
   end            end is the last chracter removed, 0 based

Note:
   sequence_name does not contain sequence comment, so be careful
   with paired end reads.
   sequence not in the trim file will be discarded
   sequence names must be in the same order as in the fastq file
```
