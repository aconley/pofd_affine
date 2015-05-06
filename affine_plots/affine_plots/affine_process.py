from __future__ import print_function, division
from collections import namedtuple
import numpy as np
from scipy.interpolate import interp1d

""" Routines related to plotting the results from
    pofd_affine fits"""


__all__ = ["read_data", "get_stats", "print_summary", "make_plots"]


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

def make_errorsnake_1D(data, nsamples=250, ninterp=100,
                       percentiles=[15.85, 84.15]):
    """ Makes 1D error snake using sampling.

    This may be a bit slow.

    This is split off from the offset, spline knots because it is 
    a log spline on f, while the others are linear."""

    if ninterp <= 2:
        raise ValueError("Ninterp must be > 2")
    if nsamples <= 2:
        raise ValueError("Nsamples must be > 2")

    logpos = np.log10(data.knots1D)
    interpvals_log = np.linspace(logpos[0], logpos[-1],
                                 ninterp, dtype=np.float32)

    n1D = len(data.knots1D)
    n1 = data.steps.shape[0]
    n2 = data.steps.shape[1]
    

    # We unfortunately have to keep track of every intermediate element
    #  so that we can later percentile them
    errorsnake_vals = np.empty((nsamples, ninterp), dtype=np.float32)

    # Recall that the knot values are already logged
    for i in range(nsamples):
        idx1 = np.random.randint(0, n1)
        idx2 = np.random.randint(0, n2)
        curr_interp = interp1d(logpos, data.steps[idx1, idx2, 0:n1D],
                               kind='cubic')
        errorsnake_vals[i, :] = curr_interp(interpvals_log)

    errorsnake_perc = np.percentile(errorsnake_vals,
                                    sorted(percentiles), axis=0)
        
    return 10.0**interpvals_log, errorsnake_perc[0, :], errorsnake_perc[1, :]


def make_color_errorsnakes(data, nsamples=150, ninterp=60,
                           percentiles=[15.85, 84.15]):
    """ Makes color param error snakes for sigma and offset.

    This may be a bit slow.

    This is split off from the 1D model because it is 
    a linear spline, while the 1D one is log on the flux density.
    Furthermore, we do these together because the conversion from
    model parameter space to physical space depends on both.
    """

    if ninterp <= 2:
        raise ValueError("Ninterp must be > 2")
    if nsamples <= 2:
        raise ValueError("Nsamples must be > 2")

    n1D = len(data.knots1D)
    nsigma = len(data.sigmaknots)
    noffset = len(data.offsetknots)
    s_interpvals = np.logspace(np.log10(data.sigmaknots[0]),
                               np.log10(data.sigmaknots[-1]),
                               ninterp, dtype=np.float32)
    o_interpvals = np.logspace(np.log10(data.offsetknots[0]),
                               np.log10(data.offsetknots[-1]),
                               ninterp, dtype=np.float32)

    n1 = data.steps.shape[0]
    n2 = data.steps.shape[1]

    # We unfortunately have to keep track of every intermediate element
    #  so that we can later percentile them
    # These will be the values in -model- space (i.e., lognorm params)
    sig_f1f2_vals = np.empty((nsamples, ninterp), dtype=np.float32)
    mu_f1f2_vals = np.empty_like(sig_f1f2_vals)

    lim1, lim2, lim3 = n1D, n1D + nsigma, n1D + nsigma + noffset
    # Part of the fun here is that the sigma and offset can be
    #  on different x scales.  So we have to do the interpolation twice
    #  and store each.
    for i in range(nsamples):
        idx1 = np.random.randint(0, n1)
        idx2 = np.random.randint(0, n2)
        s_vals = data.steps[idx1, idx2, lim1:lim2]
        s_interp = interpolator(data.sigmaknots, s_vals)
        o_vals = data.steps[idx1, idx2, lim2:lim3]
        o_interp = interpolator(data.offsetknots, o_vals)

        # Sigma f1 / f2
        curr_s = s_interp(s_interpvals)
        curr_o = o_interp(s_interpvals)
        sig_f1f2_vals[i, :] = convert_to_f1f2(curr_s, curr_o)[0]

        # Now mu
        curr_s = s_interp(o_interpvals)
        curr_o = o_interp(o_interpvals)
        mu_f1f2_vals[i, :] = convert_to_f1f2(curr_s, curr_o)[1]

    sig_f1f2_perc = np.percentile(sig_f1f2_vals, sorted(percentiles), axis=0)
    mu_f1f2_perc = np.percentile(mu_f1f2_vals, sorted(percentiles), axis=0)
    
    return s_interpvals, sig_f1f2_perc, o_interpvals, mu_f1f2_perc


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


def _other_helper(other, ax, ms, fmt, label):
    othererr = np.empty((2, other.shape[0]), dtype=np.float)
    othererr[0, :] = np.abs(other[:, 3])  # - errors
    othererr[1, :] = other[:, 2]  # + errors
    ax.errorbar(other[:, 0], other[:, 1], yerr=othererr,
                fmt=fmt, alpha=0.5, ms=ms, mew=1, fillstyle='full',
                label=label)


def band1_plot(ax, data, stats, showglenn=False, showbeth=False,
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

    # Construct central curve through points
    fit_interp = interp1d(np.log10(knots1Dplot), fit1Dplot, kind='cubic')
    plot_fvals_log = np.linspace(np.log10(knots1Dplot[0]),
                                 np.log10(knots1Dplot[-1]),
                                 npointsinterp)
    plot_fvals = 10.0**plot_fvals_log
    plot_curve = fit_interp(plot_fvals_log)

    errsnk_f, errsnk_p, errsnk_m = make_errorsnake_1D(data)

    if euclidean:
        fit1Dplot = euclideanize(knots1Dplot, fit1Dplot)
        initplot = euclideanize(knots1Dplot, initplot)
        plot_curve = euclideanize(plot_fvals, plot_curve)
        errsnk_p = euclideanize(errsnk_f, errsnk_p)
        errsnk_m = euclideanize(errsnk_f, errsnk_m)
        if showglenn:
            glenn[:, 1] = euclideanize(glenn[:, 0], glenn[:, 1])
        if showbeth:
            beth[:, 1] = euclideanize(beth[:, 0], beth[:, 1])
        if showoliver:
            oliver[:, 1] = euclideanize(oliver[:, 0], oliver[:, 1])

    # Finally, do the acutal plotting
    fiterr = np.empty((2, len(knots1Dplot)), dtype=np.float32)
    fiterr[0, :] = np.abs(unc_minus_plot)
    fiterr[1, :] = unc_plus_plot
    ax.plot(plot_fvals, plot_curve, 'r', alpha=0.3)
    ax.errorbar(knots1Dplot, fit1Dplot, yerr=fiterr, fillstyle='full',
                fmt='ro', alpha=0.7, label="This Fit")
    ax.fill_between(errsnk_f, errsnk_m, errsnk_p,
                    facecolor='grey', alpha=0.5, lw=0, edgecolor='white')

    # Oliver 2010
    if showoliver:
        _other_helper(oliver, ax, 10, 'm*', "Oliver et al. (2010)")

    # Glenn 2010
    if showglenn:
        _other_helper(glenn, ax, 5, 'bs', "Glenn et al. (2010)")

    # Bethermin 2012
    if showbeth:
        _other_helper(beth, ax, 5, 'gD', "Bethermin et al. (2012)")

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


def sigma_plot(ax, data, stats, sigma_interp, offset_interp,
               snake=None):
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

    # Add errorsnake if present
    if snake is not None:
        ax.fill_between(snake[0], snake[1][0, :], snake[1][1, :],
                        facecolor='grey', alpha=0.5, lw=0, edgecolor='white')

    
    ax.set_title("Sigma model")
    ax.set_xlabel("{0:s} Flux Density [Jy]".format(data.band1))
    ax.set_xscale('log')
    ax.set_ylabel(r"""$\sigma_{f_2 / f_1}$""")
    ax.set_xlim(0.5 * data.sigmaknots.min(),
                1.5 * data.sigmaknots.max())


def offset_plot(ax, data, stats, sigma_interp, offset_interp, snake=None):
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

    # Add errorsnake if present
    if snake is not None:
        ax.fill_between(snake[0], snake[1][0, :], snake[1][1, :],
                        facecolor='grey', alpha=0.5, lw=0, edgecolor='white')

    ax.set_title("Offset model")
    ax.set_xlabel("{0:s} Flux Density [Jy]".format(data.band1))
    ax.set_xscale('log')
    ax.set_ylabel(r"""$\left<f_2 / f_1\right>$""")
    ax.set_xlim(0.5 * data.offsetknots.min(),
                1.5 * data.offsetknots.max())


def make_plots(data, stats, showglenn=False, showbeth=False, showoliver=False,
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
    band1_plot(ax1, data, stats, showglenn=showglenn, showbeth=showbeth,
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

        # Get the errorsnake
        clr_errsnake = make_color_errorsnakes(data)
        sig_errsnake = (clr_errsnake[0], clr_errsnake[1])
        mu_errsnake = (clr_errsnake[2], clr_errsnake[3])
                                     
        if do_sigmaplot:
            ax2 = plt.subplot(gs[:, 2])
            sigma_plot(ax2, data, stats, sigma_interp, offset_interp,
                       sig_errsnake)

        if do_offsetplot:
            ax3 = plt.subplot(gs[:, 3 if do_sigmaplot else 2])
            offset_plot(ax3, data, stats, sigma_interp, offset_interp,
                        mu_errsnake)

    return f