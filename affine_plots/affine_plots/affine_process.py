from __future__ import print_function, division
from collections import namedtuple
import numpy as np
from scipy.interpolate import interp1d

""" Routines related to plotting the results from
    pofd_affine fits"""


__all__ = ["read_data", "get_stats", "print_summary", "make_plots"]


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

def make_errorsnake_1D(data, nsamples=200, ninterp=100,
                       percentiles=[15.85, 84.15]):
    """ Makes 1D error snake using sampling.

    This may be a bit slow.

    This is split off from the offset, spline knots because it is 
    a log-log spline, while the others are just log on the abcissa."""

    if ninterp <= 2:
        raise ValueError("Ninterp must be > 2")
    if nsamples <= 2:
        raise ValueError("Nsamples must be > 2")

    interpvals_log = np.linspace(np.log10(data.knots1D[0]),
                                 np.log10(data.knots1D[-1]),
                                 ninterp, dtype=np.float32)

    n1D = len(data.knots1D)
    n1 = data.steps.shape[0]
    n2 = data.steps.shape[1]
    logpos = np.log10(data.knots1D)

    # We unfortunately have to keep track of every intermediate element
    #  so that we can later percentile them
    errorsnake_vals = np.empty((nsamples, ninterp), dtype=np.float32)
    
    for i in range(nsamples):
        idx1 = np.random.randint(0, n1)
        idx2 = np.random.randint(0, n2)
        curr_interp = interp1d(logpos, np.log10(data.steps[idx1, idx2, 0:n1D]),
                               kind='cubic')
        errorsnake_vals[i, :] = curr_interp(interpvals_log)

    errorsnake_perc = np.percentile(errorsnake_vals,
                                    sorted(percentiles), axis=0)
    errorsnake_m = 10.0**errorsnake_perc[0, :]
    errorsnake_p = 10.0**errorsnake_perc[1, :]
        
    return 10.0**interpvals_log, errorsnake_p, errorsnake_m


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
               skipfirst=False, simple_errorsnake=True):
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

    # Construct the errorsnake; one way is just to plot to the
    #  uncertainties.  The better -- but more expensive way -- is
    #  to sample the curve and take the extrema
    if simple_errorsnake:
        fit_interp_plus = interp1d(np.log10(knots1Dplot),
                                   fit1Dplot + unc_plus_plot,
                                   kind='cubic')
        fit_interp_minus = interp1d(np.log10(knots1Dplot),
                                    fit1Dplot + unc_minus_plot,
                                    kind='cubic')
        errsnk_f = plot_fvals
        errsnk_p = fit_interp_plus(plot_fvals_log)
        errsnk_m = fit_interp_minus(plot_fvals_log)
    else:
        # The expensive way...
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


def sigma_plot(ax, data, stats, sigma_interp, offset_interp):
    # This will be plotted in terms of the actual sigma in
    #  f2 / f1, not the parameter.

    # TODO: Add error snake

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
    var = np.exp(2 * offset_interp(np.log10(data.sigmaknots)) +
                 sigmavals**2) * (np.exp(sigmavals**2) - 1)
    var_p = np.exp(2 * offset_interp(np.log10(data.sigmaknots)) +
                   sigmavals_p**2) * (np.exp(sigmavals_p**2) - 1)
    var_m = np.exp(2 * offset_interp(np.log10(data.sigmaknots)) +
                   sigmavals_m**2) * (np.exp(sigmavals_m**2) - 1)

    y = np.sqrt(var)
    y_minus = y - np.sqrt(var_m)
    y_plus = np.sqrt(var_p) - y

    fiterr = np.empty((2, nsigma), dtype=np.float32)
    fiterr[0, :] = y_minus
    fiterr[1, :] = y_plus
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
    offsetvals_plus = offsetvals + stats.unc_plus[idx1:idx2]
    offsetvals_minus = offsetvals + stats.unc_minus[idx1:idx2]

    # Convert from parameter to physical offset
    # Note the uncertainty in sigma is not being included, which is wrong
    # but is simpler
    mn = np.exp(offsetvals +
                0.5 * sigma_interp(np.log10(data.offsetknots))**2)
    mn_plus = np.exp(offsetvals_plus +
                     0.5 * sigma_interp(np.log10(data.offsetknots))**2)
    mn_minus = np.exp(offsetvals_minus +
                      0.5 * sigma_interp(np.log10(data.offsetknots))**2)

    y_minus = mn - mn_minus
    y_plus = mn_plus - mn

    fiterr = np.empty((2, noffset), dtype=np.float32)
    fiterr[0, :] = y_minus
    fiterr[1, :] = y_plus
    ax.plot(data.offsetknots, mn, 'r', alpha=0.3)
    ax.errorbar(data.offsetknots, mn, yerr=fiterr, fillstyle='full',
                fmt='ro', alpha=0.7, label="This Fit")
    ax.set_title("Offset model")
    ax.set_xlabel("{0:s} Flux Density [Jy]".format(data.band1))
    ax.set_xscale('log')
    ax.set_ylabel(r"""$\left<f_2 / f_1\right>$""")
    ax.set_xlim(0.5 * data.offsetknots.min(),
                1.5 * data.offsetknots.max())


def make_plots(data, stats, simple_errorsnake=True,
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
    band1_plot(ax1, data, stats, showglenn=showglenn, showbeth=showbeth,
               showoliver=showoliver, showinit=showinit,
               euclidean=euclidean, skipfirst=skipfirst,
               simple_errorsnake=simple_errorsnake)

    # The sigma and offset are plotted as constraints on f2 / f1
    #  rather than the raw parameters.  That means we need
    #  interpolation support on sigma and offset, and it makes
    #  sense to do that once and share
    if do_sigmaplot or do_offsetplot:
        n1D = len(data.knots1D)
        nsigma = len(data.sigmaknots)
        if nsigma > 1:
            sigma_interp = interp1d(np.log10(data.sigmaknots),
                                    stats.fit[n1D:(n1D+nsigma)],
                                    bounds_error=False,
                                    fill_value=stats.fit[n1D])
        else:
            # A small cheat -- just replicate the value
            sknot = np.log10(data.sigmaknots[0])
            sknotval = stats.fit[n1D]
            sigma_interp = interp1d([sknot, sknot], [sknotval, sknotval],
                                    bounds_error=False, fill_value=sknotval)
        noffset = len(data.offsetknots)
        if noffset > 1:
            offsetvals = stats.fit[(n1D+nsigma):(n1D+nsigma+noffset)]
            offset_interp = interp1d(np.log10(data.offsetknots),
                                     offsetvals, bounds_error=False,
                                     fill_value=offsetvals[0])

        else:
            oknot = np.log10(data.offsetknots[0])
            oknotval = stats.fit[n1D + nsigma]
            offset_interp = interp1d([oknot, oknot], [oknotval, oknotval],
                                     bounds_error=False,
                                     fill_value=oknotval)

        if do_sigmaplot:
            ax2 = plt.subplot(gs[:, 2])
            sigma_plot(ax2, data, stats, sigma_interp, offset_interp)

        if do_offsetplot:
            ax3 = plt.subplot(gs[:, 3 if do_sigmaplot else 2])
            offset_plot(ax3, data, stats, sigma_interp, offset_interp)

    return f
