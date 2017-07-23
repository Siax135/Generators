OBJ = pythia-6.4.28.o charges.o inelastic2.o

all : $(OBJ)
	gfortran -o inelastic2 $(OBJ)

$(OBJ) : %.o : %.f
	gfortran -c $< -o $@

clean:
	rm -f inelastic2 $(OBJ)
