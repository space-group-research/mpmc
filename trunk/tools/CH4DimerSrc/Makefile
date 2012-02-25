ch4dimer: Atom.o Molecule.o main.o
	g++ Atom.o Molecule.o main.o -o ch4dimer

main.o: main.cpp Molecule.o Atom.o
	g++ Atom.h Molecule.h main.cpp -c

Molecule.o: Molecule.h Atom.h 
	g++ Molecule.h Atom.h Molecule.cpp -c

Atom.o: Atom.h Atom.cpp
	g++ Atom.h Atom.cpp  -c
