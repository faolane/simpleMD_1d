# simple makefile for simpleMD_1d program
# Copyright (C) Fabien Brieuc - 2017
FC = gfortran
FFLAGS = -O3 -fconvert='big-endian' -Wall
FLIBS = -L/usr/local/lib -lfftw3
EXE = simpleMD_1d.exe
#OBJ = simpleMD_1d.o parameters.o thermostats.o inputOutput.o forces.o /
#energies.o ZBQ.fo
OBJ = simpleMD_1d.f90 parameters.f90 thermostats.f90 inputOutput.f90 forces.f90 energies.f90 ZBQ.f90

$(EXE): $(OBJ)
	$(FC) simpleMD_1d.f90 -o $(EXE) $(FFLAGS) $(FLIBS)

clean:
	rm -f *.mod *.exe
