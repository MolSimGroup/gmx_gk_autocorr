# preprint
This code was used for research: https://doi.org/10.1101/2023.08.28.555069.

# gmx_gk_autocorr
Calculate the autocorrelation of the pressure tensor elements from GROMACS MD simulations and the viscosity via the Green-Kubo equation. 

$\eta = \displaystyle \lim_{t\to\infty} \eta(t)  = \frac{V}{k_B T} \displaystyle \lim_{t\to\infty}  \int_0^{t} \left< P_{\alpha \beta}(t_0) \, P_{\alpha \beta}(t_0+t) \right> dt$

In addition to the off-diagonal elements, the three corresponding combinations of the diagonal pressure tensor elements ($(P_{xx}-P_{yy})/2$, $(P_{xx}-P_{zz})/2$, and $(P_{yy}-P_{zz})/2$) are used.

This code uses the output from NVT GROMACS MD simulations to calculate the autocorrelation and viscosity using the Green-Kubo equation. 

# Requirements
The code uses fast fourier transforms to compute the autocorrelation. The fast fourier transforms are implemented by using the FFTW library (https://www.fftw.org/). GROMACS is not required to run this code, but in its current version it is hard-coded to read GROMACS specific .xvg and .gro files.

# Compilation
The was compiled on a Linux machine with g++ by running e.g.

g++ -g -Wall -o main.o main.cpp readfile.cpp gkvisco.cpp -Iheaders -lfftw3 -lm

# Usage
To run the code, 2 files are necessary:
1. A .xvg file as a result from gmx energy. The order of the data in the columns is important. The order should be: "Temperature", "Pressure", "Pres-XX", "Pres-XY", "Pres-XZ","Pres-YX","Pres-YY","Pres-YZ","Pres-ZX","Pres-ZY","Pres-ZZ"
2. A .gro file representing your system size (only the last line is read to get the Volume).

To program can be executed by running the following command (change the name of the .xvg and .gro file according to your filenames (and their paths):
main.o energy.xvg conf.gro

The output files are: SELFautocorr.xvg, SELFvoltemp.xvg, and SELFvisco.xvg. These can be used for further analysis in other scripts or programs.

# Output files explained

SELFvoltemp.xvg contains the volume in nm^3 and the temperature in K.

The columns in SELFautocorr.xvg and SELVfisco.xvg list the autocorrelation or viscosity the following information:

1. time
2. ($P_{xy} + P_{yx}) / 2$
3. ($P_{xz} + P_{zx}) / 2$
4. ($P_{yz} + P_{zy}) / 2$
5. $(P_{xx}-P_{yy}) / 2$
6. $(P_{xx}-P_{zz}) / 2$
7. $(P_{yy}-P_{zz}) / 2$
8. Average of 2.-7.

The unit of 1. is in ps. The units of 2.-8. in SELFautocorr.xvg are $bar^2$. The units in 2.-8. in SELFvisco.xvg are $mPa\cdot s$.

# Example
Two example files are provided in the example_files directory. The example can be executed via "main.o exergy.xvg confout.gro".

# TO_DO
1. Add a function to also compute the viscosity curves directly.
2. Add more flexibility to read files that are not GROMACS-specific and do not need a specific order of information in the input files.
