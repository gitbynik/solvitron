all:
	gcc -o main.exe main.c consoleIN.c dsGraph.c tree.c createMat.c quantisation.c -Isrc/include -Lsrc/lib -lgsl -lgslcblas -lm