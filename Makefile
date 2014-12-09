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
	rm *.o findprimer

test:   findprimer
	findprimer -f R1.fastq.gz -p GGCACATCGATTTCTGCGAGNNNNNNNNNNNNCTCCGCTTAAGGGACT -o testR1.txt 
	findprimer -f R2.fastq.gz -p gaaaatctctagca -o testR2.txt
	join testR1.txt testR2.txt | tr " " "\t" > R1R2.txt
