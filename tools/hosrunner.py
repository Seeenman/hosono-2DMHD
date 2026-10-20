'''
Command line tool for running a hosono simulation 
'''
import argparse
from pathlib import Path
from datetime import date
from shutil import copy
from os import chdir, getcwd, environ
from subprocess import run, Popen

# create argument parser object
parser = argparse.ArgumentParser(prog='hosrunner',
                                 description='Managment script for '
                                 'setting up and running Hosono code',
                                 epilog='Written by Sean Riedel')

# add arguments to the parser object
parser.add_argument('test_problem', type=str,
                    help='Name of the test problem to run. ' \
                    'The name of the run directory will be ' \
                    'the name of the test problem with today\'s ' \
                    'date appended. ' \
                    'This will also copy the corresponding default test ' \
                    'problem parameter file from the par directory.')
parser.add_argument('--sim', type=str, nargs='+', default=[],
                    help='Simulation parameters to change. ' \
                    'Format as paramName paramValue. ' \
                    'For example, ndim 2 unsplit .true. riemannSolver \'hllc\' ')
parser.add_argument('--grid', type=str, nargs='+', default=[],
                    help='Grid parameters to change. ' \
                    'Same formatting options as simulation parameters.')
parser.add_argument('--GP', type=str, nargs='+', default=[],
                    help='Gaussian Process parameters to change. ' \
                    'Same formatting options as simulation parameters.')
parser.add_argument('--IC', type=str, nargs='+', default=[],
                    help='Initial condition parameters to change. ' \
                    'Same formatting options as simulation parameters.')
parser.add_argument('--save_to', type=str, default='./',
                    help='Optional directory to save ' \
                    'the run directory inside of relative' \
                    'to your current directory.')
parser.add_argument('--run_dir_name', type=str, 
                    help='Optional name for the run directory. '
                    'If no name given then the run directory '
                    'is automatically named based on which parameters '
                    'were changed from the default par file.')
parser.add_argument('--no_date', action='store_true',
                    help='Do not create a new directory for today\'s runs.')
parser.add_argument('--run_code', action='store_true',
                    help='Recompile and run the code.')
parser.add_argument('--time', action='store_true',
                    help='Time how long the run takes and '
                    'save the result in a file called time.txt. '
                    'Runs in the background with output saved in run.log.')
parser.add_argument('--print_output', action='store_true',
                    help='Print the output of the code to the '
                    'terminal instead of saving it in a run.log file.')
parser.add_argument('--nthreads', type=int, default=1,
                    help='Number of OpenMP threads to run the code with '
                    '(sets OMP_NUM_THREADS).')

# parse the command line arguments to the script
args = parser.parse_args()

# The parent run directory
if args.no_date:
    parent_run_dir = getcwd()+"/"+args.save_to+"/"
else:
    parent_run_dir = getcwd()+"/"+args.save_to+"/"+args.test_problem+"_"+date.today().strftime("%Y-%m-%d")+"/"

# copy the par file to the run directory 
repo_root = "~/repos/hosono-2DMHD/"
default_par_name = args.test_problem + "_default.par"
default_par_path = Path(repo_root + "par/" + default_par_name).expanduser()
if default_par_path.exists():
    # check if parent_run_dir already exists. If it doesn't then create it
    Path(parent_run_dir).mkdir(parents=True, exist_ok=True)
    copy(default_par_path, Path(parent_run_dir).expanduser())
    copied_par_file = Path(parent_run_dir+default_par_name).expanduser()
else:
    raise FileNotFoundError(f"File {default_par_path} not found. "
                     "The argument supplied to --test_problem does not "
                     "correspond to any known parameter files.")

# parse parameter file
params = {}
for raw_line in copied_par_file.read_text(encoding="utf-8").splitlines():
    line = raw_line.strip() # strip any whitespace padding

    # remove inline comments and make sure 
    # there is an equals sign present
    line = line.split("#", 1)[0].strip()
    if "=" not in line:
        continue

    key, value = line.split("=", 1)
    params[key.strip()] = value.strip()

# replace parameters
def replace_params(new_params, all_params, prefix):
    for i in range(0,len(new_params),2):
        key = prefix+'_'+new_params[i]
        val = new_params[i+1]
        if key in all_params:
            if "'" in all_params[key]:
                val = "'"+val+"'"
            elif '"' in all_params[key]:
                val = '"'+val+'"'
            all_params[key] = val
        else: 
            raise ValueError(f"Parameter \'{new_params[i]}\' not " 
                             f"found in the given parameter file with prefix \'{prefix}\'.")
    return 
replace_params(args.sim, params, 'sim')
replace_params(args.grid, params, 'grid')
replace_params(args.GP, params, 'GP')
replace_params(args.IC, params, 'IC')

# create simulation name
all_param_args = args.sim+args.grid+args.GP+args.IC
simulation_name = args.test_problem
for i in range(0,len(all_param_args),2):
    simulation_name = f"{simulation_name}-{all_param_args[i]}_{all_param_args[i+1]}"

# set the data file base name
replace_params(['dataFileBaseName','data/'+simulation_name], params, 'sim')

# write the new file
lines = Path(copied_par_file).read_text().splitlines()
new_lines = []
for raw_line in lines:
    if "=" not in raw_line:
        new_lines.append(raw_line)
        continue

    stripped = raw_line.strip()

    left, right = raw_line.split("=", 1)
    key = left.strip()

    if key not in params:
        new_lines.append(raw_line)
        continue

    # Preserve inline comment, if any
    if "#" in right:
        value_part, comment_part = right.split("#", 1)
        new_line = f"{left}= {params[key]} # {comment_part.strip()}"
    else:
        new_line = f"{left}= {params[key]}"

    new_lines.append(new_line)

if args.run_dir_name:
    run_dir_name = args.run_dir_name
else:
    run_dir_name = simulation_name

# create the subdirectory for this run and copy the new par file there
this_run_dir_path = parent_run_dir+"/"+run_dir_name+"/"
Path(this_run_dir_path).mkdir(parents=True, exist_ok=True)
new_par_file = Path(this_run_dir_path+simulation_name+".par").expanduser()
new_par_file.write_text("\n".join(new_lines) + "\n")
print(new_par_file)

# run the run directory setup shell script
chdir(this_run_dir_path)
run(Path(repo_root+"tools/run_dir_setup.sh").expanduser())

# set the number of OpenMP threads for the run
run_env = environ.copy()
if args.nthreads:
    run_env["OMP_NUM_THREADS"] = str(args.nthreads)

if args.time:
    # execute the equivalent of: (time ./hosono-2DMHD.ex run.par > run.log 2>&1) 2> time.txt &
    with open("time.txt", "w") as f:
        p = Popen(
            ["bash", "-c", 'time ./hosono-2DMHD.ex "$1" > run.log 2>&1', "bash", new_par_file],
            stderr=f,
            env=run_env,
        )
    print(p.pid)
elif args.run_code and not args.print_output:
    # execute the equivalent of the following shell command: ./hosono-2DMHD.ex run.par > run.log &
    with open("run.log", "w") as f:
        p = Popen(
            ["./hosono-2DMHD.ex", new_par_file],
            stdout=f,
            stderr=f,   # optional: also capture stderr
            env=run_env,
        )
    print(p.pid)
elif args.run_code:
    run(["./hosono-2DMHD.ex", new_par_file], env=run_env)
