# Note that := is used for immediate expansion vs deffered expansion with =
#
# Hand-rolled unit test runner, mirroring src/Makefile's style.
#
# Each test_*.f90 file in this directory is its own standalone `program` that
# links against the real src/*.f90 modules it needs (compiled straight from
# ../src, same as _build/ does for the main executable) plus assert.f90 for
# lightweight pass/fail assertions. To add a new test, drop in a new
# test_<name>.f90 next to this Makefile - it is picked up automatically.
# Add any newly-needed src module to SRC_OBJ and its module-dependency line
# below.

# Compiler
FC := gfortran

# Common compiler flags (debug-style, same as src/Makefile's FFLAGS_DEBUG).
# None of the modules under test need HDF5, so we skip those include/link flags here.
FFLAGS := -cpp -Wall -Wextra -Wimplicit-interface -Wuninitialized \
          -ffree-line-length-none -fimplicit-none -fPIC \
          -fdefault-real-8 -fdefault-double-8 \
          -ffpe-summary=invalid,zero,overflow,underflow \
          -O0 -g -ggdb -fbacktrace -fcheck=all \
          -ffpe-trap=invalid,zero,overflow \
          -finit-real=snan -finit-integer=-999999 -fmax-errors=3

SRC_DIR := ../src

# src/ modules the tests depend on
SRC_OBJ := grid.o gridBlock.o definitions.o readParamFile.o linAlgQuadPrecision.o GP.o

# local test-support modules
SUPPORT_OBJ := assert.o

# one executable per test_*.f90
TEST_SRCS := $(wildcard test_*.f90)
TEST_EXES := $(TEST_SRCS:.f90=.ex)

.PHONY: all
all: $(TEST_EXES)

.PHONY: test
test: $(TEST_EXES)
	@status=0; \
	for t in $(TEST_EXES); do \
		echo "=== running $$t ==="; \
		./$$t || status=1; \
	done; \
	exit $$status

# module dependencies (so make compiles prerequisites before dependents)
readParamFile.o: definitions.o
linAlgQuadPrecision.o: definitions.o
gridBlock.o: definitions.o
GP.o: definitions.o readParamFile.o linAlgQuadPrecision.o
grid.o: definitions.o readParamFile.o gridBlock.o
$(TEST_SRCS:.f90=.o): assert.o $(SRC_OBJ)

%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.ex: %.o $(SRC_OBJ) $(SUPPORT_OBJ)
	$(FC) $(FFLAGS) -o $@ $^

.PHONY: clean
clean:
	rm -f *.o *.mod *.ex

.PHONY: new
new: clean
	$(MAKE) test
