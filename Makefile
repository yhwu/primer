ALL: findprimer trimprimer

findprimer: primer.o functions.o
	g++ -O3 -o findprimer primer.o functions.o -lz

trimprimer: trimprimer.o functions.o
	g++ -O3 -o trimprimer trimprimer.o functions.o -lz

primer.o: primer.cpp
	g++ -c -O3 primer.cpp

trimprimer.o: trimprimer.cpp
	g++ -c -O3 trimprimer.cpp

functions.o: functions.cpp
	g++ -c -O3 functions.cpp

clean:
	rm -f *.o findprimer trimprimer
	rm -f test*.txt R1R2.txt R?.filt.fastq.gz R?cuts.txt

test:   findprimer trimprimer test.sh
	./test.sh
#	findprimer -f R1.fastq.gz -p GGCACATCGATTTCTGCGAGNNNNNNNNNNNNCTCCGCTTAAGGGACT -o testR1.txt 
#	findprimer -f R2.fastq.gz -p gaaaatctctagca -o testR2.txt
#	join testR1.txt testR2.txt | tr " " "\t" > R1R2.txt
#	awk '$5>0 && $6<5 && $10>0 && $11<4 {OFS="\t"; print $1,$4,$5,$6}' R1R2.txt > R1cuts.txt
#	awk '$5>0 && $6<5 && $10>0 && $11<4 {OFS="\t"; print $1,$9,$10,$11}' R1R2.txt > R2cuts.txt
#	trimprimer -f R1.fastq.gz -t R1cuts.txt -o R1.filt.fastq.gz
#	trimprimer -f R2.fastq.gz -t R2cuts.txt -o R2.filt.fastq.gz

