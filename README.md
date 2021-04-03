
# Setup

Clone the LpAdaptation repository

```
git clone https://github.com/breezerider/LpAdaptation_Code.git ~/LpAdaptation_Code
```

Follow these steps to build pSSAlib compatible with MATLAB. Clone the pSSAlib repository to ~/pSSAlib

```
git clone https://github.com/breezerider/pSSAlib.git ~/pSSAlib
```

and run:

```
cd ~/pSSAlib
./configure --without-libsbml --prefix=~/.local/ --disable-cli --with-pic
make clean; make -j && make install
```

# Data for Figure 1

To reproduce Figure 1, please change directory to `figure-1` and run

```
make
```

which will produce a shell script `run_fig1.sh`. Run this cripts and it will produce the true-negative and true-positive counts, that can be plotted using `plot.py`.

# Data for Figure 2

You can reproduce the results in Figure 2 by changing your working directory to `figure-2/{a,b,c}` and running 

```
make compile
```

in the respective subdirectories and running the respective scripts in the work subdirectory.
Makefiles allow one to submit individual parameter estimation runs as SLURM jobs by running:

```
make run
```

However, these MATLAB scripts can also be run interactively. First, please make sure that `LpAdaptation_Code` directory is on the MATLAB path.
Then in MATLAB just run the `.m` scripts from the respective work subdirectory. Some paramaters that you ay want to change:

* `nSampleID` -- sample identifier
* `strMethod` -- method name, 'dc' or 'cmaes'
* `strOutputPath` -- path to store the output
* `strOutputName` -- output file name prefix
* `dTimeStart` -- first time point included in the steady-state trajectory
* `dTimeStep` -- time step
* `dTimeEnd` -- last time point included in the steady-state trajectory

# Extensions

New models can be added to `models.h` and respective MATLAB driver script, following the examples presented in this repository.
