# Gaussian effective model reconstruction from contact probability matrix

The code given here implements a method to reconstruct a polymer model with harmonic pair interactions that reproduce a given contact probability matrix. There are two variants to reconstruct a Gaussian effective model:
* direct mapping;
* minimization.

## Getting Started

Follow these instructions to get a basic running example.

### Prerequisites

You will need to install [Docker](https://www.docker.com). See [the documentation](https://docs.docker.com).

### Installing

The first step is to create a docker image that contains all required compilers and applications.

```
cd dockerfiles
sudo docker build -t root/gem:16.04 16.04/.
```

Optionally you can create an additional image with user name and ID matching your local machine. Your `/home/user/` directory will be mounted at the same location within your container so that you can easily access your file system. First edit the file `dockerfiles/user/Dockerfile` and replace the user name, user ID, group name and group ID to match those on your local machine. Then:

```
cd dockerfiles
sudo docker build -t user/gem:16.04 user/.
```

In order to compile and run code of this project, you can then simply create a docker container and run the commands from inside. First create an alias:
```
alias docker-gem='sudo docker run --rm -it -w=$PWD -v $HOME:$HOME -h gem user/gem:16.04'
```

You might also consider adding directly this line to your `~/.bash_aliases` or `~/.bashrc` files. Then execute the alias to run a new container from your image in a terminal:
```
docker-gem
```

In the sequel, it is implicitly assumed that all commands are run in such a docker container.

## Minimization
Change directory to:
```
cd minimize
```

There are three implementations available:
* `minimize.cpp`.
* `minimize_thres.cpp`.
* `minimize_thres_norm.cpp`.

`minimize.cpp` requires the following arguments: `cmapfile N thres`. `cmapfile` is the path to a file containing a contact probability matrix (see examples below). `N` is the last index of the contact probability matrix to import. `thres` is the threshold to be used in the GEM mapping.

`minimize_thres.cpp` requires the following arguments: `cmapfile N thresmin thresmax dthres`. `cmapfile` is the path to a file containing a contact probability matrix. `N` is the last index of the contact probability matrix to import. `thresmin` is the minimum threshold to be used in the GEM mapping. `thresmax` is the maximum threshold. `dthres` is the threshold increment between consecutive minimizations. Only the GEM with the lowest least-square distance to the input contact probability matrix is saved.

`minimize_thres_norm.cpp` requires the following arguments: `nmapfile znorm N thresmin thresmax dthres`. `nmapfile` is the path to a file containing a contact count matrix. `znorm` is a global factor by which to divide every entry of the contact count matrix. `N` is the last index of the contact count matrix to import. `thresmin` is the minimum threshold to be used in the GEM mapping. `thresmax` is the maximum threshold. `dthres` is the threshold increment between consecutive minimizations. Only the GEM with the lowest least-square distance to the input contact probability matrix is saved.

### Compilation
The compilation process can be performed using the `minimize/bash/compile.sh` file:
```
cd minimize
bash bash/compile.sh
```
This will write an executable file named `prog`.

You may need first to edit the `compile.sh` file with a text editor to select the implementation you want to compile. You will need to comment/uncomment the lines defining the `SRC1` variable.

You may also want to modity the `KEY1` variable to choose among the following options: `gsl`, `lapack` or `mkl`. This refers to the external library used for linear algebra operations. See the source code `minimize/include/linalg.cpp` for more details.

### Running examples
#### Contact probability matrix from Brownian Dynamics simulations
Here we use a contact probability matrix computed by using configurations of a predefined GEM, sampled by Brownian Dynamics.

Go to the example directory:
```
cd examples/artificial/minimize
```

It contains the following files:
```
.
├── cmat_bd_nconf100.txt
├── cmat_th.txt
├── config.txt
├── kmat_ref.txt
├── make_matrix_plots_cmaps.yml
├── make_matrix_plots_kmaps.yml
└── run_example.sh
```

The file `kmat_ref.txt` contains the coupling matrix of the predefined GEM. The file `cmat_bd_nconf100.txt` contains the contact probability matrix obtained by sampling 100 configurations of the predefined GEM by Brownian Dynamics.

Then simply execute the script:
```
bash run_example.sh
```

The directory should now contain:
```
.
├── cmat_bd_nconf100_cmat_opt.dat
├── cmat_bd_nconf100_cmat_opt_lognorm.pdf
├── cmat_bd_nconf100.txt
├── cmat.dat
├── cmat_opt.dat
├── cmat_th.txt
├── config.txt
├── distances.dat
├── kmat.dat
├── kmat_opt.dat
├── kmat_ref_kmat_opt.dat
├── kmat_ref_kmat_opt.pdf
├── kmat_ref.txt
├── make_matrix_plots_cmaps.yml
├── make_matrix_plots_kmaps.yml
├── run_example.sh
├── sigma.dat
├── sigmainv.dat
├── sigmainv_opt.dat
└── sigma_opt.dat
```

 The contact probability matrix of the reconstructed GEM is compared to the input in `cmat_bd_nconf100_cmat_opt_lognorm.pdf`. The coupling matrix of the reconstructed GEM is compared to the original ones in the file `kmat_ref_kmat_opt.pdf`: the two matrices are very close. Note that the agreement would improve if one would use a better estimate of the predefined GEM contact probability matrix, for instance by increasing the number of Brownian Dynamics configurations used to compute the input contact probability matrix.

#### Contact probability matrix from Hi-C experiment
Here we will use a contact probability matrix coming from [published](http://dx.doi.org/10.1016/j.cell.2014.11.021) Hi-C data of human cells.

Go to the example directory:
```
cd examples/rao2014_chr8_5kbp
```

It contains the following files:
```
.
├── cmat_cmat_opt.yml
├── cmat_exp.txt
├── config.txt
├── make_matrix_plots_cmaps.yml
└── run_example.sh
```

Then simply execute the script:
```
bash run_example.sh
```

The directory should now contain:
```
.
├── cmat_cmat_opt_lognorm.pdf
├── cmat_cmat_opt.yml
├── cmat.dat
├── cmat_exp.txt
├── cmat_opt.dat
├── config.txt
├── distances.dat
├── kmat.dat
├── kmat_opt.dat
├── make_matrix_plots_cmaps.yml
├── run_example.sh
├── sigma.dat
├── sigmainv.dat
├── sigmainv_opt.dat
└── sigma_opt.dat
```

You can open `cmat_cmat_opt_lognorm.pdf` to see the comparison of the input matrix in `cmat_exp.txt` with the contact probability matrix of the reconstructed GEM in `cmat_opt.dat`.

vim: set sw=2 expandtab tabstop=2 foldcolumn=4:
vim: set spell spelllang=en_us:
