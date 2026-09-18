#!/bin/bash

test_name="test_$1"

sed "s/test_template/$test_name/g" template.f90 > $test_name.f90
touch par/$test_name.par
mkdir $test_name

cd $test_name

ln -s ../{assert.f90,Makefile,par,$test_name.f90} ./
