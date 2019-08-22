all: diploidsnp  spermsnp 

spermsnp: call_sperm_snp.cpp gzstream.cpp
	g++ -O3 -o $@  $^ -lz

diploidsnp: call_diploid_snp.cpp gzstream.cpp
	g++ -O3 -o $@  $^ -lz


clean:
	rm  diploidsnp  spermsnp
