#!/bin/bash

# save current working directory
# we do this because sometimes the run directory
# is not in the hosono git repo. 
# So to get the commit hash we need
# to make sure we are in the hosono repo.
CWD=$(pwd)
cd ~/repos/hosono/
git log --oneline --decorate | head -n 1 > $CWD/commit.txt

# compile the latest version of the code
cd ~/repos/hosono/_build_mhd/
make new
cd $CWD

# get the files we need
cp ~/repos/hosono/_build_mhd/hosono.ex ./
# cp -r ~/repos/hosono/par ./
mkdir images data
