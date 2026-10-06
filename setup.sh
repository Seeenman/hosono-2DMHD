#!/bin/bash
files=("cfl.f90" \
	"boundaryConditions.f90" \
	"choleskyDecomp.f90"     \
	"convert.f90"            \
	"definitions.f90"        \
	"eigen.f90"              \
	"eos.f90"                \
	"getFluxes.f90"          \
	"GP.f90"                 \
	"gridBlock.f90"          \
	"grid.f90"               \
	"gridState.f90"          \
	"initialCondition.f90"   \
	"linAlgQuadPrecision.f90"\
	"mhd_driver.f90"         \
	"output.f90"             \
	"readParamFile.f90"      \
	"reconstruct.f90"        \
	"riemannSolver.f90"      \
	"simulation.f90"         \
	"timeStep.f90"           \
	"Makefile"               \
)

build_dir=$1
build_dir_is_new=false
[[ ! -d $build_dir ]] && mkdir $build_dir && build_dir_is_new=true

if [ $build_dir_is_new ]; then
    for f in "${files[@]}"; do
	ln -s "$(pwd)"/src/"$f" "$build_dir/"
    done
fi
