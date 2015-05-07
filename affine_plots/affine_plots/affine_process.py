from __future__ import print_function, division
from collections import namedtuple
import numpy as np
from scipy.interpolate import interp1d

""" Routines related to plotting the results from
    pofd_affine fits"""


__all__ = ["read_data", "get_stats", "print_summary", "make_plots",
           "read_errorsnake"]


class interpolator:
    """ Wrapper around interp1d that handles edges correctly"""

    # An issue here is that none of the scipy splines actually
    # match what we are using (GSL cspline).  However, it is clear
    # from experimentation that the cubic spline from scipy is a lot
    # more prone to wild excursions than the GSL one, so we stick
    # to that
    def __init__(self, knotpos, knotval):
        self._npos = len(knotpos)
        self._pos = knotpos.copy()
        self._val = knotval.copy()

        if self._npos == 1:
            self._interp = None
            if not np.isscalar(self._val):
                self._val = self._val[0]
        elif self._npos == 2:
            self._interp = interp1d(self._pos, self._val, kind='slinear')
        elif self._npos > 2:
            self._interp = interp1d(self._pos, self._val, kind='quadratic')
        else:
            raise ValueError("Not enough points")

    def __call__(self, xvals):
        if self._interp is None:
            retvals = self._val * np.ones_like(xvals)
        else:
            retvals = np.empty_like(xvals)
            leftpos = self._pos[0]
            leftval = self._val[0]
            rightpos = self._pos[-1]
            rightval = self._val[-1]
            for idx, x in enumerate(xvals):
                if x < leftpos:
                    retvals[idx] = leftval
                elif x > rightpos:
                    retvals[idx] = rightval
                else:
                    retvals[idx] = self._interp(x)

        return retvals   

def convert_to_f1f2(sigma, offset):
    """ Convert from lognormal param space to f1 / f2"""

    var_f1of2 = np.exp(2 * offset + sigma**2) * np.expm1(sigma**2)
    mu_f1of2 = np.exp(offset + 0.5 * sigma**2)
    return np.sqrt(var_f1of2), mu_f1of2
    

def find_bandname(h5filename):
    """ Tries to guess the band from h5file name.

    Basically, search the filename for things that look
    like P[SML]W, so this is specific to Herschel/SPIRE."""

    import re
    regex = re.compile("P[SML]W", re.IGNORECASE)
    return [m.upper() for m in regex.findall(h5filename)][:2]


def read_data(h5file):
    """ Read data from HDF5 file"""
    import os.path
    import h5py

    if not os.path.isfile(h5file):
        raise IOError("Can't open HDF5 input file {0:s}".format(h5file))

    is1D = True
    with h5py.File(h5file, 'r') as f:
        steps = f['Chains/Chain'][...]
        like = f['Chains/Likelihood'][...]
        knots1D = f['LikelihoodParams/Model/KnotPositions'][...]
        initvals = f['ParamInfo/InitialPosition'][...]
        param_names = [b.decode() for b in
                       f['ParamInfo'].attrs['ParamNames']]
        if 'LikelihoodParams/Model/SigmaKnotPositions' in f:
            is1D = False
            sigmaknots = f['LikelihoodParams/Model/SigmaKnotPositions'][...]
            offsetknots = f['LikelihoodParams/Model/OffsetKnotPositions'][...]

    if is1D:
        mcmc_data1D = namedtuple('mcmc_data1D',
                                 ['steps', 'like', 'knots1D', 'allknots',
                                  'initvals', 'param_names', 'band',
                                  'npars', 'nparsfit', 'nbonus'])
        bands = find_bandname(h5file)
        band = bands[0] if len(bands) > 0 else "Unknown"
        data = mcmc_data1D(steps, like, knots1D, knots1D, initvals,
                           param_names, band, steps.shape[2],
                           steps.shape[2] - 2, 2)
    else:
        mcmc_data2D = namedtuple('mcmc_data2D',
                                 ['steps', 'like', 'knots1D', 'sigmaknots',
                                  'offsetknots', 'allknots',
                                  'initvals', 'param_names',
                                  'band1', 'band2', 'npars',
                                  'nparsfit', 'nbonus'])
        bands = find_bandname(h5file)
        band1 = bands[0] if len(bands) >= 2 else "Unknown"
        band2 = bands[1] if len(bands) >= 2 else "Unknown"
        allknots = np.hstack([knots1D, sigmaknots, offsetknots])
        data = mcmc_data2D(steps, like, knots1D, sigmaknots, offsetknots,
                           allknots, initvals, param_names, band1, band2,
                           steps.shape[2], steps.shape[2] - 4, 4)
    return data


def read_errorsnake(h5file):
    """ Read errorsnake from HDF5 file"""
    import os.path
    import h5py

    if not os.path.isfile(h5file):
        raise IOError("Can't open HDF5 input file {0:s}".format(h5file))

    is1D = True
    with h5py.File(h5file, 'r') as f:
        flux1D = f['Flux'][...]
        mean1D = f['MeanLog10Counts'][...]
        median1D = f['MedianLog10Counts'][...]
        snake1D_0p683 = f['Log10CountsSnakep=0.683'][...]
        snake1D_0p954 = f['Log10CountsSnakep=0.954'][...]
        if 'FluxSigma' in f:
            is1D = False

    if is1D:
        snake_data1D = namedtuple('snake_data1D',
                                  ['flux1D', 'mean1D', 'median1D',
                                   'snake1D_0p683', 'snake1D_0p954'])
        data = snake_data1D(flux1D, mean1D, median1D, snake1D_0p683,
                            snake1D_0p954)
    else:
        raise NotImplementedError("2D error snake not yet implemented")

    return data


def get_stats(data, thin=5, percentiles=[15.85, 84.15]):
    """ Compute statistics on mcmc results"""

    if thin <= 0:
        raise ValueError("Invalid thinning value: {0:d}".format(thin))

    bidx = np.unravel_index(data.like.argmax(), data.like.shape)

    stats_type = namedtuple('stats', ['fit', 'unc_plus', 'unc_minus',
                                      'best', 'bestlike'])

    stats = stats_type(np.empty(data.npars, dtype=np.float32),
                       np.empty(data.npars, dtype=np.float32),
                       np.empty(data.npars, dtype=np.float32),
                       data.steps[bidx[0], bidx[1], :], data.like[bidx])

    for i in range(data.npars):
        vals = data.steps[:, ::thin, i]
        stats.fit[i] = vals.mean()
        perc = np.percentile(vals, percentiles)
        stats.unc_plus[i] = perc[1] - stats.fit[i]
        stats.unc_minus[i] = perc[0] - stats.fit[i]

    return stats


def print_summary(data, stats):
    """ Prints summary statistics"""

    # Print
    print("#### Model ####")
    print("%-11s  %-6s  %5s %6s %6s %6s %6s" %
          ("#Param", "Knot", "Init", "Best", "Fit", "Unc+", "Unc-"))
    for i in range(data.nparsfit - 2):  # -2 for sigma multipliers
        print("%-11s  %6g  %5.2f %6.3f %6.3f %+6.3f %+6.3f" %
              (data.param_names[i], data.allknots[i], data.initvals[i],
               stats.best[i], stats.fit[i],
               stats.unc_plus[i], stats.unc_minus[i]))

    print("#### Sigma multipliers ####")
    print("%-11s  %5s %6s %6s %6s %6s" %
          ("#Param", "Init", "Best", "Fit", "Unc+", "Unc-"))
    for i in range(data.nparsfit - 2, data.nparsfit):
        print("%-11s  %5.2f %6.3f %6.3f %+6.3f %+6.3f" %
              (data.param_names[i], data.initvals[i],
               stats.best[i], stats.fit[i],
               stats.unc_plus[i], stats.unc_minus[i]))

    print("#### Bonus params ####")
    print("%-11s  %8s %8s %7s %7s" %
          ("#Param", "Best", "Fit", "Unc+", "Unc-"))
    for i in range(data.nparsfit, data.npars):
        print("%-11s  %8.3f %8.3f %+7.3f %+7.3f" %
              (data.param_names[i], stats.best[i], stats.fit[i],
               stats.unc_plus[i], stats.unc_minus[i]))

    print("Best log like: %0.2f" % stats.bestlike)


def euclideanize(knotpos, knotvals):
    return knotvals + 2.5 * np.log10(knotpos)


def otherdata_helper(other, ax, ms, fmt, label):
    othererr = np.empty((2, other.shape[0]), dtype=np.float)
    othererr[0, :] = np.abs(other[:, 3])  # - errors
    othererr[1, :] = other[:, 2]  # + errors
    ax.errorbar(other[:, 0], other[:, 1], yerr=othererr,
                fmt=fmt, alpha=0.5, ms=ms, mew=1, fillstyle='full',
                label=label)


def band1_plot(ax, data, stats, errorsnake=None,
               showglenn=False, showbeth=False,
               showoliver=False, showinit=False, euclidean=False,
               skipfirst=False):
    """ Does plot for 1D band model"""
    from pkg_resources import resource_filename

    npointsinterp = 100

    # Read in comparison data if needed
    bandname = data.band1 if hasattr(data, 'band1') else data.band
    if showglenn:
        glennfile = resource_filename(__name__,
                                      'resources/power_nocfirb_{0:s}.txt')
        glenn = np.loadtxt(glennfile.format(bandname.lower()))
    if showbeth:
        bethfile = resource_filename(__name__,
                                     'resources/bethermin2012_{0:s}.txt')
        beth = np.loadtxt(bethfile.format(bandname.lower()))
    if showoliver:
        oliverfile = resource_filename(__name__,
                                       'resources/oliver2010_{0:s}.txt')
        oliver = np.loadtxt(oliverfile.format(bandname.lower()))

    # Peel off points we will actually plot
    if skipfirst:
        sidx = 1
    else:
        sidx = 0
    n1D = len(data.knots1D)
    knots1Dplot = data.knots1D[sidx:].copy()
    fit1Dplot = stats.fit[sidx:n1D].copy()
    unc_plus_plot = stats.unc_plus[sidx:n1D].copy()
    unc_minus_plot = stats.unc_minus[sidx:n1D].copy()
    bestplot = stats.best[sidx:n1D].copy()
    initplot = data.initvals[sidx:n1D].copy()

    # Error snake if present
    if errorsnake is not None:
        errsnake_f = errorsnake.flux1D.copy()
        errsnake_med = errorsnake.median1D.copy()
        errsnake_mn = errorsnake.mean1D.copy()
        errsnake_m1 = errorsnake.snake1D_0p683[:, 0].copy()
        errsnake_p1 = errorsnake.snake1D_0p683[:, 1].copy()
        errsnake_m2 = errorsnake.snake1D_0p954[:, 0].copy()
        errsnake_p2 = errorsnake.snake1D_0p954[:, 1].copy()

    if euclidean:
        fit1Dplot = euclideanize(knots1Dplot, fit1Dplot)
        initplot = euclideanize(knots1Dplot, initplot)
        if showglenn:
            glenn[:, 1] = euclideanize(glenn[:, 0], glenn[:, 1])
        if showbeth:
            beth[:, 1] = euclideanize(beth[:, 0], beth[:, 1])
        if showoliver:
            oliver[:, 1] = euclideanize(oliver[:, 0], oliver[:, 1])
        if errorsnake is not None:
            errsnake_mn = euclideanize(errsnake_f, errsnake_mn)
            errsnake_med = euclideanize(errsnake_f, errsnake_med)
            errsnake_m1 = euclideanize(errsnake_f, errsnake_m1)
            errsnake_p1 = euclideanize(errsnake_f, errsnake_p1)
            errsnake_m2 = euclideanize(errsnake_f, errsnake_m2)
            errsnake_p2 = euclideanize(errsnake_f, errsnake_p2)

    # Finally, do the acutal plotting
    fiterr = np.empty((2, len(knots1Dplot)), dtype=np.float32)
    fiterr[0, :] = np.abs(unc_minus_plot)
    fiterr[1, :] = unc_plus_plot
    ax.errorbar(knots1Dplot, fit1Dplot, yerr=fiterr, fillstyle='full',
                fmt='ro', alpha=0.7, label="This Fit")

    if errorsnake is not None:
        ax.plot(errsnake_f, errsnake_mn, 'r', alpha=0.3)
        ax.fill_between(errsnake_f, errsnake_m2, errsnake_p2,
                        facecolor='grey', alpha=0.2, lw=0, edgecolor='grey')
        ax.fill_between(errsnake_f, errsnake_m1, errsnake_p1,
                        facecolor='grey', alpha=0.4, lw=0, edgecolor='grey')

    # Oliver 2010
    if showoliver:
        otherdata_helper(oliver, ax, 10, 'm*', "Oliver et al. (2010)")

    # Glenn 2010
    if showglenn:
        otherdata_helper(glenn, ax, 5, 'bs', "Glenn et al. (2010)")

    # Bethermin 2012
    if showbeth:
        otherdata_helper(beth, ax, 5, 'gD', "Bethermin et al. (2012)")

    # Plot initial values
    if showinit:
        ax.plot(knots, initplot, 'ro', ms=1.5, alpha=0.5,
                label="Initial values")

    ax.set_title("{0:s} model".format(bandname))
    ax.set_xlabel('Flux Density [Jy]')
    if euclidean:
        ylab = r"""$\log_{10}\,S^{\,2.5}dN/dS$  [Jy$^{1.5}$ deg$^{-2}$]"""
    else:
        ylab = r"""$\log_{10}\,dN/dS$  [Jy$^{-1}$ deg$^{-2}$]"""
    ax.set_ylabel(ylab)

    ax.set_xscale('log')
    ax.set_xlim(0.8 * knots1Dplot.min(), 1.2 * knots1Dplot.max())
    if euclidean:
        ax.legend(loc=3, fontsize='small')
    else:
        ax.legend(fontsize='small')


def sigma_plot(ax, data, stats, sigma_interp, offset_interp):
    # This will be plotted in terms of the actual sigma in
    #  f2 / f1, not the parameter.

    n1D = len(data.knots1D)
    nsigma = len(data.sigmaknots)
    if nsigma < 2:
        raise ValueError("Don't call sigma_plot with one sigma knot")

    sigmavals = stats.fit[n1D:(n1D + nsigma)]
    sigmavals_p = sigmavals + stats.unc_plus[n1D:(n1D + nsigma)]
    sigmavals_m = sigmavals + stats.unc_minus[n1D:(n1D + nsigma)]

    # Convert from parameter to physical sigma
    # Note the uncertainty in mu is not being included, which is wrong
    # but is simpler
    y = convert_to_f1f2(sigmavals, offset_interp(data.sigmaknots))[0]
    y_m = convert_to_f1f2(sigmavals_m, offset_interp(data.sigmaknots))[0]
    y_p = convert_to_f1f2(sigmavals_p, offset_interp(data.sigmaknots))[0]

    fiterr = np.empty((2, nsigma), dtype=np.float32)
    fiterr[0, :] = y - y_m
    fiterr[1, :] = y_p - y
    ax.plot(data.sigmaknots, y, 'r', alpha=0.3)
    ax.errorbar(data.sigmaknots, y, yerr=fiterr, fillstyle='full',
                fmt='ro', alpha=0.7, label="This Fit")

    
    ax.set_title("Sigma model")
    ax.set_xlabel("{0:s} Flux Density [Jy]".format(data.band1))
    ax.set_xscale('log')
    ax.set_ylabel(r"""$\sigma_{f_2 / f_1}$""")
    ax.set_xlim(0.5 * data.sigmaknots.min(),
                1.5 * data.sigmaknots.max())


def offset_plot(ax, data, stats, sigma_interp, offset_interp):
    # This will be plotted in terms of the actual mean of f2/f1
    # rather than the model mu parameter.

    n1D = len(data.knots1D)
    nsigma = len(data.sigmaknots)
    noffset = len(data.offsetknots)
    if noffset < 2:
        raise ValueError("Don't call offset_plot with one offset knot")

    idx1 = n1D + nsigma
    idx2 = idx1 + noffset
    offsetvals = stats.fit[idx1:idx2]
    offsetvals_p = offsetvals + stats.unc_plus[idx1:idx2]
    offsetvals_m = offsetvals + stats.unc_minus[idx1:idx2]

    # Convert from parameter to physical offset
    # Note the uncertainty in sigma is not being included, which is wrong
    # but is simpler
    y = convert_to_f1f2(sigma_interp(data.offsetknots), offsetvals)[1]
    y_m = convert_to_f1f2(sigma_interp(data.offsetknots), offsetvals_p)[1]
    y_p = convert_to_f1f2(sigma_interp(data.offsetknots), offsetvals_m)[1]

    fiterr = np.empty((2, noffset), dtype=np.float32)
    fiterr[0, :] = y - y_m
    fiterr[1, :] = y_p - y
    ax.plot(data.offsetknots, y, 'r', alpha=0.3)
    ax.errorbar(data.offsetknots, y, yerr=fiterr, fillstyle='full',
                fmt='ro', alpha=0.7, label="This Fit")

    ax.set_title("Offset model")
    ax.set_xlabel("{0:s} Flux Density [Jy]".format(data.band1))
    ax.set_xscale('log')
    ax.set_ylabel(r"""$\left<f_2 / f_1\right>$""")
    ax.set_xlim(0.5 * data.offsetknots.min(),
                1.5 * data.offsetknots.max())


def make_plots(data, stats, errorsnake=None,
               showglenn=False, showbeth=False, showoliver=False,
               showinit=False, euclidean=False, skipfirst=False):
    """ Plot results of a fit"""

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    # Use the same machinery for 1D and 2D, even though
    #  it's ridiculous overkill in 1D.
    # We are potentially plotting up to 3 panels in the 2D case: 
    #  1D fit, sigma, offset
    # However, sigma and offset could be single points, in which
    #  case plotting is pointless
    do_sigmaplot = False
    do_offsetplot = False
    is2D = hasattr(data, 'sigmaknots')
    if is2D:
        do_sigmaplot = len(data.sigmaknots) > 1
        do_offsetplot = len(data.offsetknots) > 1

    # The band 1 uses 2 panels, others use at most 1
    n_panels = 2 + (1 if do_sigmaplot else 0) + (1 if do_offsetplot else 0)
    figsize = {2: (10, 5), 3: (13, 5), 4: (16, 5)}

    f = plt.figure(figsize=figsize[n_panels])
    gs = gridspec.GridSpec(1, n_panels)
    gs.update(left=0.08, right=0.92, wspace=0.35)

    if is2D:
        # 1D doesn't really need a supertitle
        f.suptitle("{0:s} vs. {1:s}".format(data.band1, data.band2))

    ax1 = plt.subplot(gs[:, 0:2])
    band1_plot(ax1, data, stats, errorsnake=errorsnake,
               showglenn=showglenn, showbeth=showbeth,
               showoliver=showoliver, showinit=showinit,
               euclidean=euclidean, skipfirst=skipfirst)

    # The sigma and offset are plotted as constraints on f2 / f1
    #  rather than the raw parameters.  That means we need
    #  interpolation support on sigma and offset, and it makes
    #  sense to do that once and share
    if do_sigmaplot or do_offsetplot:
        n1D = len(data.knots1D)
        nsigma = len(data.sigmaknots)
        sigma_interp = interpolator(data.sigmaknots,
                                    stats.fit[n1D:(n1D+nsigma)])

        noffset = len(data.offsetknots)
        lim1 = n1D + nsigma
        lim2 = lim1 + noffset
        offset_interp = interpolator(data.offsetknots,
                                     stats.fit[lim1:lim2])

                                     
        if do_sigmaplot:
            ax2 = plt.subplot(gs[:, 2])
            sigma_plot(ax2, data, stats, sigma_interp, offset_interp)

        if do_offsetplot:
            ax3 = plt.subplot(gs[:, 3 if do_sigmaplot else 2])
            offset_plot(ax3, data, stats, sigma_interp, offset_interp)

    return f
