#!/usr/bin/env python3

""" Print the results of a 2D P(D) fit"""

from collections import namedtuple
import os.path
import h5py
import numpy
from scipy.interpolate import interp1d

def find_bandname(h5filename):
    """ Determine band from h5file name.  Assumed to be of the form
    band1band2_something.h5"""

    spl = h5filename.split('_')
    if len(spl) < 2:
        errstr = "Don't know how to get band from {0:s}"
        raise ValueError(errstr.format(h5filename))
    if len(spl[0]) != 6:
        errstr = "Don't know how to get 2 bands from {0:s}"
        raise ValueError(errstr.format(spl[0]))
    bands = {'PSW', 'PMW', 'PLW'}
    bandname1 = spl[0][0:3].upper()
    bandname2 = spl[0][3:].upper()
    if bandname1 not in bands:
        raise ValueError("Unknown band {0:s}".format(bandname1))
    if bandname2 not in bands:
        raise ValueError("Unknown band {0:s}".format(bandname2))
    return bandname1, bandname2

def read_data(h5file):

    if not os.path.isfile(h5file):
        raise IOError("Can't open HDF5 input file {0:s}".format(h5file))

    mcmc_data = namedtuple('mcmc_data', ['steps', 'like', 'knots1D',
                                         'sigmaknots', 'offsetknots',
                                         'initvals', 'param_names',
                                         'band1', 'band2'])

    with h5py.File(h5file, 'r') as f: 
        steps = f['Chains/Chain'][...]
        like = f['Chains/Likelihood'][...]
        knots1D = f['LikelihoodParams/Model/KnotPositions'][...]
        sigmaknots = f['LikelihoodParams/Model/SigmaKnotPositions'][...]
        offsetknots = f['LikelihoodParams/Model/OffsetKnotPositions'][...]
        initvals =  f['ParamInfo/InitialPosition'][...]
        param_names = [b.decode() for b in f['ParamInfo'].attrs['ParamNames']]

    band1, band2 = find_bandname(h5file)
    return mcmc_data(steps, like, knots1D, sigmaknots, offsetknots,
                     initvals, param_names, band1, band2)


def get_stats(data, thin=5):

    if thin <= 0:
        raise ValueError("Invalid thinning value: {0:d}".format(thin))

    bidx = numpy.unravel_index(data.like.argmax(), data.like.shape)
    
    stats_type = namedtuple('stats', ['fit', 'unc_plus', 'unc_minus',
                                      'best', 'bestlike'])

    npars = data.steps.shape[2]
    stats = stats_type(numpy.empty(npars, dtype=numpy.float32),
                       numpy.empty(npars, dtype=numpy.float32),
                       numpy.empty(npars, dtype=numpy.float32),
                       data.steps[bidx[0], bidx[1], :], data.like[bidx])
                       
    for i in range(npars):
        vals = data.steps[:, ::thin, i]
        stats.fit[i] = vals.mean()
        perc = numpy.percentile(vals, [15.85, 84.15])  # 1 sigma
        stats.unc_plus[i] = perc[1] - stats.fit[i]
        stats.unc_minus[i] = perc[0] - stats.fit[i]

    return stats

def print_summary(data, stats):
    all_knots = numpy.hstack((data.knots1D, data.sigmaknots, data.offsetknots))
    nknots = len(data.knots1D)
    nsigmaknots = len(data.sigmaknots)
    noffsetknots = len(data.offsetknots)
    npars = data.steps.shape[2]
    nparfit = data.steps.shape[2] - 4 # -4 for bonus params

    # Print
    print("#### Model ####")
    print("%-11s  %-6s  %5s %6s %6s %6s %6s" %\
          ("#Param", "Knot","Init","Best","Fit","Unc+", "Unc-"))
    for i in range(nparfit - 2):  # -2 for sigma multipliers
        print("%-11s  %6g  %5.2f %6.3f %6.3f %+6.3f %+6.3f" %\
              (data.param_names[i], all_knots[i], data.initvals[i],
               stats.best[i], stats.fit[i],
               stats.unc_plus[i], stats.unc_minus[i]))

    print("#### Sigma multipliers ####")
    print("%-11s  %5s %6s %6s %6s %6s" %\
          ("#Param", "Init", "Best", "Fit", "Unc+", "Unc-"))
    for i in range(nparfit - 2, nparfit):
        print("%-11s  %5.2f %6.3f %6.3f %+6.3f %+6.3f" %\
              (data.param_names[i], data.initvals[i],
               stats.best[i], stats.fit[i],
               stats.unc_plus[i], stats.unc_minus[i]))

    print("#### Bonus params ####")
    print("%-11s  %8s %8s %7s %7s" %\
          ("#Param", "Best", "Fit", "Unc+", "Unc-"))
    for i in range(nparfit, npars):
        print("%-11s  %8.3f %8.3f %+7.3f %+7.3f" %\
              (data.param_names[i], stats.best[i], stats.fit[i],
               stats.unc_plus[i], stats.unc_minus[i]))
        
    print("Best log like: %0.2f" % stats.bestlike)

    
def band1_plot(ax, data, stats, showglenn=False, showbeth=False,
               showoliver=False, showinit=False, euclidean=False,
               skipfirst=False):

    npointsinterp = 100
    
    n1D = len(data.knots1D)

    # Read in comparison data if needed
    otherdir = '/Users/aconley/data/PofD/fits/1D/other/'
    if showglenn:
        glennfile = os.path.join(otherdir, 'power_nocfirb_{0:s}.txt')
        glenn = numpy.loadtxt(glennfile.format(data.band1.lower()))
    if showbeth:
        bethfile = os.path.join(otherdir, 'bethermin2012_{0:s}.txt')
        beth = numpy.loadtxt(bethfile.format(data.band1.lower()))
    if showoliver:
        oliverfile = os.path.join(otherdir, 'oliver2010_{0:s}.txt')
        oliver = numpy.loadtxt(oliverfile.format(data.band1.lower()))
        
    
    # Peel off points we will actually plot
    if skipfirst:
        sidx = 1
    else:
        sidx = 0
    knots1Dplot = data.knots1D[sidx:].copy()
    fit1Dplot = stats.fit[sidx:n1D].copy()
    unc_plus_plot = stats.unc_plus[sidx:n1D].copy()
    unc_minus_plot = stats.unc_minus[sidx:n1D].copy()
    bestplot = stats.best[sidx:n1D].copy()
    initplot = data.initvals[sidx:n1D].copy()
    
    # Construct error snake using log-log interpolation.
    # The error snake is not really being done correctly here,
    # since that would be expensive.
    fit_interp = interp1d(numpy.log10(knots1Dplot), fit1Dplot,
                          kind='cubic')
    fit_interp_plus = interp1d(numpy.log10(knots1Dplot), 
                               fit1Dplot + unc_plus_plot,
                               kind='cubic')
    fit_interp_minus = interp1d(numpy.log10(knots1Dplot), 
                                fit1Dplot + unc_minus_plot,
                                kind='cubic')
    plot_fvals_log = numpy.linspace(numpy.log10(knots1Dplot[0]), 
                                    numpy.log10(knots1Dplot[-1]),
                                    npointsinterp)
    plot_curve = fit_interp(plot_fvals_log)
    plot_curve_plus = fit_interp_plus(plot_fvals_log)
    plot_curve_minus = fit_interp_minus(plot_fvals_log)
    plot_fvals = 10**plot_fvals_log

    if euclidean:
        fit1Dplot += 2.5 * numpy.log10(knots1Dplot)
        initplot += 2.5 * numpy.log10(knots1Dplot)
        plot_curve += 2.5 * plot_fvals_log
        plot_curve_plus += 2.5 * plot_fvals_log
        plot_curve_minus += 2.5 * plot_fvals_log
        if showglenn:
            glenn[:, 1] += 2.5 * numpy.log10(glenn[:, 0])
        if showbeth:
            beth[:, 1] += 2.5 * numpy.log10(beth[:, 0])
        if showoliver:
            oliver[:, 1] += 2.5 * numpy.log10(oliver[:, 0])

    # Finally, do the acutal plotting
    fiterr = numpy.empty((2, len(knots1Dplot)), dtype=numpy.float32)
    fiterr[0, :] = numpy.abs(unc_minus_plot)
    fiterr[1, :] = unc_plus_plot
    ax.plot(plot_fvals, plot_curve, 'r', alpha=0.3)
    ax.errorbar(knots1Dplot, fit1Dplot, yerr=fiterr, fillstyle='full',
                fmt='ro', alpha=0.7, label="This Fit")
    ax.fill_between(plot_fvals, plot_curve_minus, plot_curve_plus,
                    facecolor='grey', alpha=0.5, lw=0, edgecolor='white')

    # Oliver 2010
    if showoliver:
        olivererr = numpy.empty((2, oliver.shape[0]), dtype=numpy.float)
        olivererr[0, :] = numpy.abs(oliver[:, 3])  # - errors
        olivererr[1, :] = oliver[:, 2]  # + errors
        ax.errorbar(oliver[:, 0], oliver[:, 1], yerr=olivererr,
                    fmt='m*', alpha=0.5, ms=10, mew=1, fillstyle='full',
                    label="Oliver et al. (2010)")
        
    # Glenn 2010
    if showglenn:
        glennerr = numpy.empty((2, glenn.shape[0]), dtype=numpy.float)
        glennerr[0, :] = numpy.abs(glenn[:, 3])  # - errors
        glennerr[1, :] = glenn[:, 2]  # + errors
        ax.errorbar(glenn[:, 0], glenn[:, 1], yerr=glennerr,
                    fmt='bs', alpha=0.5, mew=1, fillstyle='full',
                    label="Glenn et al. (2010)")
            
    # Bethermin 2012
    if showbeth:
        betherr = numpy.empty((2, beth.shape[0]), dtype=numpy.float)
        betherr[0, :] = numpy.abs(beth[:, 3])  # - errors
        betherr[1, :] = beth[:, 2]  # + errors
        ax.errorbar(beth[:, 0], beth[:, 1], yerr=betherr,
                    fmt='gD', alpha=0.5, mew=1, 
                    fillstyle='full', label="Bethermin et al. (2012)")

    # Plot initial values
    if showinit:
        ax.plot(knots, initplot, 'ro', ms=1.5, alpha=0.5,
                label="Initial values")


    ax.set_title("{0:s} model".format(data.band1))
    ax.set_xlabel('Flux Density [Jy]')
    if euclidean:
        ax.set_ylabel(r"""$\log_{10}\,S^{\,2.5}dN/dS$  [Jy$^{1.5}$ deg$^{-2}$]""")
    else:
        ax.set_ylabel(r"""$\log_{10}\,dN/dS$  [Jy$^{-1}$ deg$^{-2}$]""")

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
    sigmavals_plus = sigmavals + stats.unc_plus[n1D:(n1D + nsigma)]
    sigmavals_minus = sigmavals + stats.unc_minus[n1D:(n1D + nsigma)]
    
    # Convert from parameter to physical sigma
    # Note the uncertainty in mu is not being included, which is wrong
    # but is simpler
    var = numpy.exp(2 * offset_interp(numpy.log10(data.sigmaknots)) +
                    sigmavals**2) * (numpy.exp(sigmavals**2) - 1)
    var_plus = numpy.exp(2 * offset_interp(numpy.log10(data.sigmaknots)) +
                         sigmavals_plus**2) * \
                         (numpy.exp(sigmavals_plus**2) - 1)
    var_minus = numpy.exp(2 * offset_interp(numpy.log10(data.sigmaknots)) +
                          sigmavals_minus**2) * \
                          (numpy.exp(sigmavals_minus**2) - 1)

    y = numpy.sqrt(var)
    y_minus = y - numpy.sqrt(var_minus)
    y_plus = numpy.sqrt(var_plus) - y
    
    fiterr = numpy.empty((2, nsigma), dtype=numpy.float32)
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

    # TODO: Add error snake

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
    mn = numpy.exp(offsetvals +
                   0.5 * sigma_interp(numpy.log10(data.offsetknots))**2)
    mn_plus = numpy.exp(offsetvals_plus +
                        0.5 * sigma_interp(numpy.log10(data.offsetknots))**2)
    mn_minus = numpy.exp(offsetvals_minus +
                         0.5 * sigma_interp(numpy.log10(data.offsetknots))**2)
    
    y_minus = mn - mn_minus
    y_plus = mn_plus - mn
    
    fiterr = numpy.empty((2, noffset), dtype=numpy.float32)
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
    
            
def plot_results(data, stats, showglenn=False, showbeth=False,
                 showoliver=False, showinit=False, euclidean=False,
                 skipfirst=False):
    """ Plot results of fit"""

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    
    # In the general case, we plot 3 panels:
    #  1D fit, sigma, offset
    # But sigma and offset could be single points, in which
    #  case plotting is pointless
    do_sigmaplot = len(data.sigmaknots) > 1
    do_offsetplot = len(data.offsetknots) > 1

    # The band 1 uses 2 panels
    n_panels = 2 + (1 if do_sigmaplot else 0) + (1 if do_offsetplot else 0)
    figsize = {2: (10,5), 3: (13,5), 4: (16,5)}

    f = plt.figure(figsize=figsize[n_panels])
    gs = gridspec.GridSpec(1, n_panels)
    gs.update(left=0.08, right=0.92, wspace=0.35)
    
    f.suptitle("{0:s} vs. {1:s}".format(data.band1, data.band2))

    ax1 = plt.subplot(gs[:, 0:2])
    band1_plot(ax1, data, stats, showglenn=showglenn, showbeth=showbeth,
               showoliver=showoliver, showinit=showinit, euclidean=euclidean,
               skipfirst=skipfirst)

    # The sigma and offset are plotted as constraints on f2 / f1
    #  rather than the raw parameters.  That means we need
    #  interpolation support on sigma and offset, and it makes
    #  sense to do that once and share
    if do_sigmaplot or do_offsetplot:
        n1D = len(data.knots1D)
        nsigma = len(data.sigmaknots)
        if nsigma > 1:
            sigma_interp = interp1d(numpy.log10(data.sigmaknots),
                                    stats.fit[n1D:(n1D+nsigma)],
                                    bounds_error=False,
                                    fill_value=stats.fit[n1D])
        else:
            # A small cheat -- just replicate the value
            sknot = numpy.log10(data.sigmaknots[0])
            sknotval = stats.fit[n1D]
            sigma_interp = interp1d([sknot, sknot], [sknotval, sknotval],
                                    bounds_error=False, fill_value=sknotval)
        noffset = len(data.offsetknots)
        if noffset > 1:
            offsetvals = stats.fit[(n1D+nsigma):(n1D+nsigma+noffset)]
            offset_interp = interp1d(numpy.log10(data.offsetknots),
                                     offsetvals, bounds_error=False,
                                     fill_value=offsetvals[0])
                                     
        else:
            oknot = numpy.log10(data.offsetknots[0])
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
    
if __name__ == "__main__":
    import argparse

    desc = """Plot and print the results of a P(D, D) fit"""
    
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('h5file', action='store', 
                        help="HDF5 file to load")
    parser.add_argument('--thin', action='store', type=int,
                        default=5, help="Thinning")
    parser.add_argument('--beth', action='store_true', default=False,
                        help="Show Bethermin 2012 data")
    parser.add_argument('--euclidean', action='store_true', default=False,
                        help="Plot as euclidean normalized counts")
    parser.add_argument('--glenn', action='store_true', default=False,
                        help="Show Glenn 2010 data")
    parser.add_argument('--noplot', action='store_true', default=False,
                        help="Don't make plots, just print statistics")
    parser.add_argument('--oliver', action='store_true', default=False,
                        help="Show Oliver 2010 data")
    parser.add_argument('-o', '--outfile', action='store', default=None,
                        help="Output PDF file")
    parser.add_argument('--plotinit', action='store_true', default=False,
                        help="Plot initial values")
    parser.add_argument('--skipfirst', action='store_true', default=False,
                        help="Don't plot first point")
    
    parse_results = parser.parse_args()  # Runs on sys.argv by default

    data = read_data(parse_results.h5file)
    stats = get_stats(data, thin=parse_results.thin)
    print_summary(data, stats)
    f = plot_results(data, stats, showglenn=parse_results.glenn,
                     showbeth=parse_results.beth,
                     showoliver=parse_results.oliver,
                     showinit=parse_results.plotinit,
                     euclidean=parse_results.euclidean,
                     skipfirst=parse_results.skipfirst)

    if parse_results.outfile is None:
        f.savefig(os.path.splitext(parse_results.h5file)[0] + '.pdf')
    else:
        f.savefig(outfile)
