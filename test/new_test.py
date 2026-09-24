#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("name")
args = parser.parse_args()

test_name = f"test_{args.name}"

template = Path("template.f90").read_text()
Path(f"{test_name}.f90").write_text(template.replace("test_template", test_name))
Path(f"par/{test_name}.par").touch()
os.mkdir(test_name)

os.chdir(test_name)

for f in ["assert.f90", "Makefile", "par", f"{test_name}.f90"]:
    os.symlink(f"../{f}", f"./{f}")
