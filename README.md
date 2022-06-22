# DIABLO3T
## Reinventing the triply-periodic version of DIABLO with modern Fortran and the 2DECOMP&FFT library

**This code is at the early stages of development and is not yet ready for production use!**
### What is done so far
- Time-stepper following original diablo_per approach (RK3 + Crank-Nicolson)
- Output of scalar fields to HDF5 restart files
- Initial condition prescription routines
- Reading of old-style diablo_per `input.dat` file

### Still to be implemented
- Output out of flow field
- Reading of output files to use as input
- Calculation of statistics (and output of these)
- Buoyancy coupling
- Background gradients
- Custom forcing

## Quick start
This version of DIABLO relies on the 2DECOMP&FFT library, which needs installing before DIABLO3T can be compiled.
Simply clone the library into your home directory by:

``
git clone https://github.com/numericalalgorithmsgroup/2decomp_fft.git
``

Edit some variables in `src/Makefile.inc` to build the 2DECOMP library with your favourite FFT library.
I recommend setting the flags

``
OPTIONS=-DDOUBLE_PREC -DEVEN
FFT=fftw3
``

The `-DEVEN` flag will mean that the grid must divide exactly by the number of processes in each dimension, but should be slightly more optimised.
You will also need to set `FFTW_PATH` when using `FFT=fftw3`.
For a standard Ubuntu installation using apt, this will be `FFTW_PATH=/usr`, and on a cluster you can check the relevant FFT module information to find the correct path.

Run `make` in the base directory of the 2DECOMP library to build it, then you should be able to build DIABLO3T using the Makefile.