#Makefile for Taylor with GCC-10 compilers

SHELL = /bin/bash

DEBUG = false#{true, false}  (NOTE: no whitespace at end)

FC = gfortran

#*************************************************

ifeq ($(DEBUG), true)
FCFLAGS = -c -O0 -Wall -fcheck=all -fbacktrace -std=f2018 \
		  -fno-unsafe-math-optimizations \
		  -frounding-math -fsignaling-nans
else
FCFLAGS = -c -O3 -std=f2008ts \
		  -funroll-loops \
		  -flto=auto     \
		  -fno-unsafe-math-optimizations \
		  -frounding-math -fsignaling-nans
endif

SRC_DIR = ../src/core ../src
BUILD_DIR = .
INSTALL_DIR = ../bin
#EXECUTABLE = ligan

vpath %.f90 $(SRC_DIR)
#vpath $(EXECUTABLE) $(INSTALL_DIR)

INCS = #-I$(SRC_DIR) 
#LIBS = 

OBJ_NAMES = \
		  constants_m   \
		  strings_m     \
		  utils_math_m  \
		  qsort_m       \
		  vector_m      \
		  M_CLI2        \
		  trajectory_reader_m  \
		  simbox_m      \
		  specs_m       \
		  atmcfg_m      \
		  config_io_m

OBJECTS = $(addsuffix .o,$(OBJ_NAMES))

.PHONY : all
all : bin2txt

.PHONY : bin2txt
bin2txt : $(INSTALL_DIR)/bin2txt

$(INSTALL_DIR)/bin2txt : $(OBJECTS) bin2txt.o
	if [ -d "$(INSTALL_DIR)" ]; then mkdir -p "$(INSTALL_DIR)"; fi
	$(FC) $(INCS) $^ $(LIBS) -o $@

%.o : %.f90
	$(FC) $(INCS) $(FCFLAGS) $^ -o $@

.PHONY : clean clobber
clean:
	-rm -rf *.o *.mod
clobber:
	-rm -rf $(INSTALL_DIR)/bin2txt *.o *.mod
