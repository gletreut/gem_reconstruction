import os,sys,shutil
import numpy as np
import yaml

############################################################################
##  PARAMETERS
############################################################################
origin = os.path.dirname(os.path.relpath(sys.argv[0]))
name = os.path.basename(sys.argv[0])
name = os.path.splitext(name)[0]

sys.path.append(os.path.join(origin,"py"))
from functions import make_distance_plot

############################################################################
## MAIN
############################################################################
# arguments
if not (len(sys.argv) > 2):
    print "syntax: <param file> <file1> <file2> ... <fileN>"
    sys.exit()
paramfile = os.path.relpath(sys.argv[1])
if not (os.path.isfile(paramfile)):
    print "File {} does not exist!".format(paramfile)
    sys.exit()
print "{:<20s}{:<s}".format("paramfile",paramfile)
if not (os.path.splitext(paramfile)[1] == ".yml"):
    raise ValueError("param file must be a yaml file!")
fin = open(paramfile,"r")
params = yaml.load(fin)
fin.close()

# loop in the target directories
tfiles=sys.argv[2:]
for tfile in tfiles:
    if not (os.path.isfile(tfile)):
        print "file {} does not exist. Skipping...".format(tfile)
        continue

    # copy param file for records
    tdir = os.path.dirname(tfile)
    filename = os.path.basename(paramfile)
    dest = os.path.join(tdir,filename)
    if (os.path.realpath(dest) != os.path.realpath(paramfile)):
        shutil.copyfile(paramfile,dest)

    fileout = make_distance_plot(tfile, **params['make_distance_plot_args'])
    print "{:<20s}{:<s}".format("fileout",fileout)

