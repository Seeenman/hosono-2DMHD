#!/bin/bash

repo_root="$HOME/repos/hosono-2DMHD/"
build_dir="_build"

# save current working directory
# we do this because sometimes the run directory
# is not in the hosono git repo. 
# So to get the commit hash we need
# to make sure we are in the hosono repo.
CWD=$(pwd)
cd $repo_root
git log --oneline --decorate | head -n 1 > $CWD/commit.txt

# compile the latest version of the code
cd $repo_root/$build_dir
make clean && make release
cd $CWD

# get the files we need
cp $repo_root/$build_dir/hosono-2DMHD.ex ./
mkdir images data
