import os,sys,shutil
import numpy as np
import yaml

############################################################################
##  PARAMETERS
############################################################################
origin = os.path.dirname(os.path.relpath(sys.argv[0]))
name = os.path.basename(sys.argv[0])
name = os.path.splitext(name)[0]

#paramfile = os.path.join(origin,"roles","{}.yml".format(name))
#if not (os.path.isfile(paramfile)):
#    print "File {} does not exist!".format(paramfile)
#    sys.exit()
#print "{:<20s}{:<s}".format("paramfile",paramfile)
#
#fin = open(paramfile,"r")
#params = yaml.load(fin)
#fin.close()

sys.path.append(os.path.join(origin,"py"))
from functions import read_map, make_matrix_plot, make_matrix_plot_two, make_matrix_plot_multi, make_matrix_compare_multi

############################################################################
## MAIN
############################################################################
# arguments
if not (len(sys.argv) > 2):
    print "syntax: <param file> <tdir1> <tdir2> ... <tdirN>"
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

tdirs=sys.argv[2:]
for tdir in tdirs:
    if (not os.path.isdir(tdir)):
        print "directory does not exist: {}".format(tdir)
        sys.exit()
    print "{:<20s}{:<s}".format("tdir",tdir)

    # copy param file for records
    filename = os.path.basename(paramfile)
    dest = os.path.join(tdir,filename)
    try:
        shutil.copyfile(paramfile,dest)
    except shutil.Error:
        pass

# loop in the target directories
for tdir in tdirs:
    # make single plots
    for basename in params['single_matrices']:
        f = os.path.join(tdir,basename)
        if not (os.path.isfile(f)):
            print "file {} does not exist. Skipping...".format(f)
            continue
        fileout = make_matrix_plot(f,**params['make_matrix_plot_arg'])
        print "{:<20s}{:<s}".format("fileout",fileout)

    # make multi plots
    if (len(params['multi_matrices']) == 1):
        file_exp = params['multi_matrices'][0][0]
        file_pred = params['multi_matrices'][0][1]
        fs=[]
        basenames = [file_exp,file_pred]
        for basename in basenames:
            f = os.path.join(tdir,basename)
            if not (os.path.isfile(f)):
                sys.exit("some files in the list do not exist: {}".format(" ".join(basenames)))
            fs.append(f)
        fileout = make_matrix_plot_two(fs[0],fs[1],**params['make_matrix_plot_two_arg'])
        print "{:<20s}{:<s}".format("fileout",fileout)
        fileout = make_matrix_compare_multi(fs, **params['make_matrix_compare_arg'])
        print "{:<20s}{:<s}".format("fileout",fileout)
    else:
        for basenames in params['multi_matrices']:
            skip=False
            fs = []
            for basename in basenames:
                f = os.path.join(tdir,basename)
                if not (os.path.isfile(f)):
                    skip=True
                    print "some files in the list do not exist: {}".format(" ".join(basenames))
                fs.append(f)
            if (skip):
                continue
            fileout = make_matrix_plot_multi(fs,**params['make_matrix_plot_multi_arg'])
            print "{:<20s}{:<s}".format("fileout",fileout)
            fileout = make_matrix_compare_multi(fs, **params['make_matrix_compare_arg'])
            print "{:<20s}{:<s}".format("fileout",fileout)

