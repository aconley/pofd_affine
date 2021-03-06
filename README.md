pofd_affine
===========

This is a package to measure number counts from confused
images using the P(D) formalism and an affine-invariant
MCMC.  Both one band (i.e., number counts) and 2D
(number counts with color information) models are supported.
It was developed with sub-mm/mm astronomy in mind, but
is in principle useful for other problems.  However, it does
assume sub-mm/mm conventions, such as per-beam normalization
of the input data.

### Installation

Installation is via the standard UNIX `configure` and
`make`. pofd_affine depends on the following packages:
* [FFTW](http://www.fftw.org/).  Generating a 
  [FFTW wisdom](http://www.fftw.org/fftw-wisdom.1.html)
  file will speed up the code signficantly.
* [cfitsio](http://heasarc.gsfc.nasa.gov/fitsio/)
* [hdf5](http://www.hdfgroup.org/HDF5/)
* The [GNU scientific library](http://www.gnu.org/software/gsl/)
* MPI. There are many versions of MPI, but 
   [OpenMPI](http://www.open-mpi.org/) is what I used. 
* If you want to use the built in tests, you'll need 
   [googletest](http://code.google.com/p/googletest/)
It may be necessary to tell configure where to look for these
libraries -- see `configure --help`.

Note that this package requires C++11 support from your compiler;
it may be necessary to add compiler flags to turn this on.

In addition, pofd_affine comes with an optional Python script 
to print summaries and make plots of the results.  This
can be installed by going into the affine_plots subdirectory
and executing `python setup.py install`.

This requires the Python libraries:
* [h5py](http://www.h5py.org/)
* [matplotlib](http://matplotlib.org/)
* [numpy](http://www.numpy.org/)
* [scipy](http://www.scipy.org)
Although probably the best way to install them is to use `pip`,
if for some reason they aren't already present.

### Documentation

There isn't much.  Doxygen style documentation can be
generated using

	make docs

but that isn't really all that helpful.   All of the command
line routines have some built in documentation using `--help`:

	pofd_affine_getPD --help

Note that some routines -- including `pofd_affine_mcmc` -- have to
be run using MPI.

There are some built in tests (unit and otherwise) that can be built
if you have used --enable-test when configuring, and then do

        make check

There are also some examples, which can be built if --enable-examples
is set while configuring.

The basic idea is to create a specfile that tells the code about
the data, an initfile that specifies the model and initial position,
and then run something like

      mpiexec pofd_affine_mcmc -H initfile specfile outfile.h5

Note that this can take a while to run.  Summary statistics are
then created using

    pofd_affine_getStats outfile.h5 stats.h5

and can be plotted if the optional python package is installed using

    make_affine_plots stats.h5  

All of these routines have numerous additional options.

### Branches

The master branch contains the P(D) code.  The base branch has
only code related to the affine-invariant MCMC.  

### References
* The original P(D) paper is [Scheuer (1957)](http://dx.doi.org/10.1017/S0305004100032825),
* A more recent and easier to read discussion is
  [Patanchon et al. (2009)](http://dx.doi.org/10.1088/0004-637X/707/2/1750),
  which gives results from BLAST.
* The most recent P(D) results from SPIRE can be found at
  [Glenn et al. (2010)](http://dx.doi.org/10.1111/j.1365-2966.2010.17781.x).
* The affine invariant MCMC method used in this code is detailed
  in [Goodman and Weare (2010)](http://msp.berkeley.edu/camcos/2010/5-1/camcos-v5-n1-p04-p.pdf),
  and the method of parallelizing it in
  [Foreman-Macket et al. (2012)](http://arxiv.org/abs/1202.3665); the MCMC
  routines here are to first order a re-implementation of the
  [emcee](https://github.com/dfm/emcee) package in c++ with lots of
  P(D) routines added.
