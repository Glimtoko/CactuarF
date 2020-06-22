# Fortran compiler and base flags
FC := gfortran
FFLAGS := -O3

# Object list
all_objects =  \
	riemann_sampler.o \
	riemann_solvers.o \
	riemann.o \
	flux.o \
	cactaur.o

# Installation directory
INSTALL_DIR := /prod/CactaurF/bin/

# Expected library locations
PROD := /usr/local/prod
TIO_DIR := $(PROD)/typhonio
SILO_DIR := $(PROD)/silo-4.10.2-bsd
HDF5_DIR := $(PROD)/hdf5_oldapi

# Include and link directories
INCLUDES :=
LINKS :=
RPATH :=
PREPRO :=
NEED_HDF = NO


# Check for libraries
TIO_CHECK := $(TIO_DIR)/include/tio_access.h

$(info Checking for TyphonIO in $(TIO_DIR))
ifneq ("$(wildcard $(TIO_CHECK))", "")
  $(info  TIO Found)
  INCLUDES += -I$(TIO_DIR)/include
  LINKS += -L$(TIO_DIR)/lib -ltyphonio_f -ltyphonio
  PREPRO += -DTIO

  NEED_HDF = YES
else
  $(info TIO Not found)
endif


SILO_CHECK := $(SILO_DIR)/include/silo_f9x.inc

$(info )
$(info Checking for SILO in $(SILO_DIR))
ifneq ("$(wildcard $(SILO_CHECK))", "")
  $(info  SILO Found)
  INCLUDES += -I$(SILO_DIR)/include
  LINKS += -L$(SILO_DIR)/lib -lsiloh5
  PREPRO += -DSILO

  NEED_HDF = YES
else
  $(info SILO Not found)
endif


ifeq ($(NEED_HDF), YES)
  $(info )
  $(info HDF5 Library required)
  HDF5_CHECK := $(HDF5_DIR)/include/h5a.mod

  $(info Checking for HDF5 in $(HDF5_DIR))
  ifneq ("$(wildcard $(HDF5_CHECK))", "")
    $(info  HDF5 Found)
    INCLUDES += -I$(HDF5_DIR)/include
    LINKS += -L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5
    RPATH += -Wl,-rpath $(PROD)/hdf5/lib
  endif
endif

$(info )
$(info CactaurF Environment::)
$(info Fortran compiler: $(FC))
$(info Base compiler flags: $(FFLAGS))
$(info Preprocessor flags:$(PREPRO))
$(info )

.SUFFIXES:
.SUFFIXES: .o .f95 .F95
.PHONY: clean install

.f95.o:
	$(FC) $(FFLAGS) $(INCLUDES) $< -c

.F95.o:
	$(FC) $(FFLAGS) $(PREPRO) $(INCLUDES) $< -c

cactuarf: $(all_objects)
	$(FC) $(FFLAGS) $(all_objects) $(LINKS) $(RPATH) -o cactuarf

clean:
	-rm cactuarf
	-rm *.o
	-rm *.mod

install: cactuarf
	@echo Installing to $(INSTALL_DIR)
	cp cactuarf $(INSTALL_DIR)
	chmod +x $(INSTALL_DIR)/cactuarf


# Dependencies
prog.o: data.o memory_manipulation.o mesh_reader.o sod.o lagstep_fi.o text_output.o

sod.o: data.o	# For input and file unit

text_output.o: write_tio.o write_silo.o

