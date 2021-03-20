FC = gfortran
.SUFFIXES : .f90

FFLAGS = -g -Wall -fbounds-check -Wuninitialized -ffpe-trap=invalid,zero,overflow

TARGET = ../bin/GPSCALPOS
OBJECTS = main.o mod_variable.o time_util.o nav_msg.o  satellite_position.o print_list.o compute_solution.o
MOD_FILE = mod_variable.mod time_util.mod nav_msg.mod  satellite_position.mod print_list.mod compute_solution.mod

%.o: %.f90
	$(COMPILE.f) $(OUTPUT_OPTION) $<

%.mod: %.f90 %.o
	@:

${TARGET} : ${OBJECTS}
	$(FC) -o $@ ${OBJECTS} $(FFLAGS)

main.o : ${MOD_FILE}


.PHONY: clean run

clean:
	rm -f ${OBJECTS} ${MOD_FILE}

run:
	../bin/GPSCALPOS