test: matrix.f90 compute_solution.f90 main.f90
	gfortran -g -Wall -fbounds-check -O -Wuninitialized -ffpe-trap=invalid,zero,overflow -o test compute_solution.f90 main.f90
