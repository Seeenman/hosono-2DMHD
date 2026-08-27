# Note that := is used for immediate expansion vs deffered expansion with =

# Compiler
FC := gfortran

# Target
TARGET := hosono-2DMHD.ex

# ==============================================================
# Since we are using hdf5, we need to specify the paths to the 
# include files and library files for hdf5 (sometimes have to 
# do this even if using the h5fc hdf5 fortran compiler wrapper
# for some reason).
# ==============================================================
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    HDF5_INCLUDE_PATH := /usr/include
    HDF5_LIBRARY_PATH := /usr/lib
else ifeq ($(UNAME_S),Darwin)
    HDF5_INCLUDE_PATH := /opt/homebrew/include
    HDF5_LIBRARY_PATH := /opt/homebrew/lib
endif

# C Preprocessor (CPP) and include flags
CPPFLAGS := -I$(HDF5_INCLUDE_PATH)

# Linker flags
LDFLAGS  := -L$(HDF5_LIBRARY_PATH)

# libraries to link against
LDLIBS   := -lhdf5_fortran -lhdf5

# Common compiler flags
FFLAGS_COMMON := -cpp -Wall -Wextra -Wimplicit-interface -Wuninitialized \
                 -ffree-line-length-none -fimplicit-none -fPIC \
                 -fdefault-real-8 -fdefault-double-8 -fopenmp \
                 -ffpe-summary=invalid,zero,overflow,underflow

# Debug flags
FFLAGS_DEBUG := -O0 -g -ggdb -fbacktrace -fcheck=all \
                -ffpe-trap=invalid,zero,overflow \
                -finit-real=snan -finit-integer=-999999

# Optimization flags
FFLAGS_RELEASE := -O3 -march=native -funroll-loops

# Default build mode (?= means assign only if not already set)
MODE ?= debug

ifeq ($(MODE),debug)
    FFLAGS := $(FFLAGS_COMMON) $(FFLAGS_DEBUG)
else ifeq ($(MODE),release)
    FFLAGS := $(FFLAGS_COMMON) $(FFLAGS_RELEASE)
else
    $(error MODE must be 'debug' or 'release')
endif

OBJ := definitions.o readParamFile.o simulation.o grid.o \
       boundaryConditions.o eos.o convert.o initialCondition.o output.o \
       mhd_driver.o

.PHONY: all
all: $(TARGET)

.PHONY: debug
debug:
	$(MAKE) MODE=debug

.PHONY: release
release:
	$(MAKE) MODE=release

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)
	@echo "============="
	@echo "===SUCCESS==="
	@echo "============="

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(TARGET) *.o *.mod

.PHONY: tidy
tidy:
	rm -f *.o *.mod

.PHONY: new # equivalent to running make clean && make
new: clean
	$(MAKE)
