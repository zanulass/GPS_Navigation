test: matrix.f90 compute_solution.f90 main.f90
	gfortran -o test matrix.f90 compute_solution.f90 main.f90
