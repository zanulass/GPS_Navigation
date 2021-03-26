FC = gfortran
.SUFFIXES : .f90

FFLAGS = -g -Wall -fbounds-check -Wuninitialized -ffpe-trap=invalid,zero,overflow

TARGET = ../bin/GPSCALPOS

OBJECTS = main.o \
	WGS84_constants.o \
	exec_conditions.o \
	compute_solution.o \
	time_util.o \
	coordinate_util.o \
  navigation_message.o \
	print_util.o \
	satellite_position.o \
	correction_parameters.o \
	compute_receiver_position.o





MOD_FILE = WGS84_constants.mod \
	exec_conditions.mod \
	compute_solution.mod \
	time_util.mod \
	coordinate_util.mod \
	navigation_message.mod \
	print_util.mod \
	satellite_position.mod \
	correction_parameters.mod \
	compute_receiver_position.mod





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
	rm -f ../tmp/*
	../bin/GPSCALPOS