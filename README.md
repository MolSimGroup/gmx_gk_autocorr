# gmx_gk_autocorr
Calculate the autocorrelation of the pressure tensor elements from GROMACS MD simulations, which is needed to compute the viscosity via the Green-Kubo equation. 

$\eta = \displaystyle \lim_{t\to\infty} \eta(t)  = \frac{V}{k_B T} \displaystyle \lim_{t\to\infty}  \int_0^{t} \left< P_{\alpha \beta}(t_0) \, P_{\alpha \beta}(t_0+t) \right> dt$

In addition to the offdiagonal elements, the three corresponding combinations of the diagonal pressure tensor elements ($(P_{xx}-P_{yy})/2$, $(P_{xx}-P_{zz})/2$, and $(P_{yy}-P_{zz})/2$) are used.

This code uses the output from NVT GROMACS MD simulations to calculate the autocorrelation using the Green-Kubo equation. 

# Requirements
The code uses fast fourier transforms to compute the autocorrelation. The fast fourier transforms are implemented by using the FFTW library (https://www.fftw.org/). GROMACS is not required to run this code, but in its current version it is hard coded to read GROMACS specific .xvg and .gro files.

# Compilation
The was compiled on a linux machine with g++ by running e.g.

g++ -g -Wall -o main.o main.cpp readfile.cpp gkvisco.cpp -Iheaders -lfftw3 -lm

# Usage
To run the code, 2 files are necessary:
1. A .xvg file as a result from gmx energy. The order of the data in the columns is important. The order should be: @ s0 legend "Temperature", "Pressure", "Pres-XX", "Pres-XY", "Pres-XZ","Pres-YX","Pres-YY","Pres-YZ","Pres-ZX","Pres-ZY","Pres-ZZ"
2. A .gro file representing your system size (only the last line is read to get the Volume).

To program can be executed by running the following command (change the name of the .xvg and .gro file according to your filenames (and their paths):
main.o energy.xvg conf.gro

The output files are: SELFautocorr.xvg and SELFvoltemp.xvg These can be used for further analysis in other scripts or programs. (Currently, also a SELFvisco.xvg is written out, not including the diagonal combinations. This still has to be expanded).

# TO_DO
1. Add a function to also compute the viscosity curves directly.
2. Add more flexibility to read files that are not GROMACS-specific and do not need a specific order of information in the input files.
