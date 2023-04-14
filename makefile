CC = g++
FLAG = -Wall -Wextra -pedantic -O3 -pthread -march=native -fopenmp
LINKER = -lntl -lfplll -lgmp

all:Generator.exe test.exe

#val:Generator.exe test.exe
#	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose ./Generator.exe -P_order 10

Generator.exe:Generator.o Base.o MNSGenerator.o Norm.o FindM.o FindRoots.o Misc.o Bound.o Meta.o Poly.o
	$(CC) $^ $(FLAG) -o $@ $(LINKER)

test.exe:test.o Base.o MNSGenerator.o Norm.o FindM.o FindRoots.o Misc.o Bound.o Meta.o Poly.o
	$(CC) $^ $(FLAG) -o $@ $(LINKER)

Generator.o:Generator.cpp 
	$(CC) $< $(FLAG) -c -o $@ $(LINKER)

test.o:test.cpp 
	$(CC) $< $(FLAG) -c -o $@ $(LINKER)

MNSGenerator.o:MNSGenerator.cpp MNSGenerator.h
	$(CC) $< $(FLAG) -c -o $@ $(LINKER)

Norm.o:Norm.cpp Norm.h
	$(CC) $< $(FLAG) -c -o $@ $(LINKER)

Base.o:Base.cpp Base.h
	$(CC) $< $(FLAG) -c -o $@ $(LINKER)

FindM.o:FindM.cpp FindM.h
	$(CC) $< $(FLAG) -c -o $@ $(LINKER)

FindRoots.o:FindRoots.cpp FindRoots.h
	$(CC) $< $(FLAG) -c -o $@ $(LINKER)

Misc.o:Misc.cpp Misc.h
	$(CC) $< $(FLAG) -c -o $@ $(LINKER)

Bound.o:Bound.cpp Bound.h
	$(CC) $< $(FLAG) -c -o $@ $(LINKER)

Meta.o:Meta.cpp Meta.h
	$(CC) $< $(FLAG) -c -o $@ $(LINKER)

Poly.o:Poly.cpp Poly.h
	$(CC) $< $(FLAG) -c -o $@ $(LINKER)