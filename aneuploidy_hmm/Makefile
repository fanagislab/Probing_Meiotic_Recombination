all: cnvhmm

cnvhmm: hmm.cpp main.cpp
	g++ -O3 -static -o $@  $^ 

clean:
	rm cnvhmm

install:
	cp cnvhmm /gpfsdata/Analysis/fanwei/exe
