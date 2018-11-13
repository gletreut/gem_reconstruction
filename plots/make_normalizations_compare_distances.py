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
from functions import make_distance_superimposition_plot

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

# Detect normalization for each file
tfiles=sys.argv[2:]
dist_dict = {}
parse_key = params['parse_key']
for tfile in tfiles:
    ## existence
    if not (os.path.isfile(tfile)):
        print "file {} does not exist. Skipping...".format(tfile)
        continue

    ## get value of normalization used
    tdir = os.path.dirname(os.path.realpath(tfile))
    name = os.path.basename(tdir)
    tab = name.split('_')
    val = None
    for el in tab:
        ind = el.find(parse_key)
        if (ind > -1):
            val = float(el[ind+len(parse_key):])
            break

    if (val != None):
        dist_dict[val] = tfile

# select norms to show
selection = params['selection']
if (selection == None):
    selection = dist_dict.keys()

distfiles=[]
norms=[]
for s in selection:
    if not (s in dist_dict.keys()):
        continue
    norms.append(s)
    distfiles.append(dist_dict[s])

if (len(distfiles) == 0):
    raise ValueError("distance files list is empty!")


# plot
nij_max = params['nij_max']
tdir = os.path.dirname(os.path.dirname(os.path.relpath(tfiles[0])))
fileout = make_distance_superimposition_plot(norms=norms, distfiles=distfiles, nij_max=nij_max, tdir=tdir, **params['make_distance_superimposition_plot_args'])
print "{:<20s}{:<s}".format("fileout",fileout)

