from __future__ import print_function, division
from collections import namedtuple
import numpy as np
from scipy.interpolate import interp1d

""" Routines related to plotting the results from
    pofd_affine fits"""


__all__ = ["read_data", "print_summary", "make_plots"]


def find_bandname(h5filename):
    """ Tries to guess the band from h5file name.

    Basically, search the filename for things that look
    like P[SML]W, so this is specific to Herschel/SPIRE."""

    import re
    regex = re.compile("P[SML]W", re.IGNORECASE)
    return [m.upper() for m in regex.findall(h5filename)][:2]


def read_data(h5file):
    """ Read data from HDF5 stats file"""
    import os.path
    import h5py

    if not os.path.isfile(h5file):
        raise IOError("Can't open HDF5 input file {0:s}".format(h5file))

    with h5py.File(h5file, 'r') as f:
        modelType = f['ModelType'][...][0].decode()
        if f['NDims'][0] == 1:
            stats1D = namedtuple('stats1D',
                                 ['Band', 'ModelType', 'ParamNames',
                                  'KnotPositions', 'ParamsInit',
                                  'BestLikelihood', 'ParamsBest', 
                                  'ParamsMean', 'Params1Sig', 
                                  'Snake1DFlux', 'Snake1DMean', 
                                  'Snake1D1Sig', 'Snake1D2Sig'])

            # Values of actual params
            param_names = [p.decode() for p in f['Params/ParamNames']]
            knot_positions = f['Params/KnotPositions'][...]
            best_like = f['Params/BestLikelihood'][0]
            params_init = f['Params/ParamsInit'][...]
            params_best = f['Params/ParamsBest'][...]
            params_mean = f['Params/ParamsMean'][...]
            params_1sig = f['Params/Paramsp=0.683'][...]
            params_1sig -= params_mean[:, np.newaxis]

            bands = find_bandname(h5file)
            band = bands[0] if len(bands) > 0 else "Unknown"

            # Smooth plotting curves (errorsnakes, etc.)
            curve_f = f['Curves/FluxDensity'][...]
            curve_mean = f['Curves/Log10CountsMean'][...]
            curve_1sig = f['Curves/Log10Countsp=0.683'][...]
            curve_2sig = f['Curves/Log10Countsp=0.954'][...]
        
            data = stats1D(band, modelType, param_names, knot_positions,
                           params_init, best_like, params_best, 
                           params_mean, params_1sig, curve_f, 
                           curve_mean, curve_1sig, curve_2sig)
        elif f['NDims'][0] == 2:
            stats2D = namedtuple('stats2D',
                                 ['Band1', 'Band2', 'ModelType',
                                  'AreRatios', 'ParamNames',
                                  'KnotPositions', 'SigmaKnotPositions',
                                  'OffsetKnotPositions', 'ParamsInit',
                                  'BestLikelihood', 'ParamsBest', 
                                  'ParamsMean', 'Params1Sig', 
                                  'Snake1DFlux', 'Snake1DMean', 
                                  'Snake1D1Sig', 'Snake1D2Sig',
                                  'SnakeSigmaFlux', 'SnakeSigmaMean',
                                  'SnakeSigma1Sig', 'SnakeSigma2Sig',
                                  'SnakeOffsetFlux','SnakeOffsetMean', 
                                  'SnakeOffset1Sig', 'SnakeOffset2Sig'])

            # Values of actual params
            are_ratios = f['AreF2F1Ratios'][0]
            param_names = [p.decode() for p in f['Params/ParamNames']]
            knot_positions = f['Params/KnotPositions'][...]
            sigma_positions = f['Params/SigmaKnotPositions'][...]
            offset_positions = f['Params/OffsetKnotPositions'][...]
            params_init = f['Params/ParamsInit'][...]
            best_like = f['Params/BestLikelihood'][0]
            params_best = f['Params/ParamsBest'][...]
            params_mean = f['Params/ParamsMean'][...]
            params_1sig = f['Params/Paramsp=0.683'][...]
            params_1sig -= params_mean[:, np.newaxis]

            # Smooth plotting curves (errorsnakes, etc.)
            curve1D_f = f['Curves/Band1/FluxDensity'][...]
            curve1D_mean = f['Curves/Band1/Log10CountsMean'][...]
            curve1D_1sig = f['Curves/Band1/Log10Countsp=0.683'][...]
            curve1D_2sig = f['Curves/Band1/Log10Countsp=0.954'][...]
            curveSigma_f = f['Curves/Color/SigmaFluxDensity'][...]
            curveSigma_mean = f['Curves/Color/SigmaMean'][...]
            curveSigma_1sig = f['Curves/Color/Sigmap=0.683'][...]
            curveSigma_2sig = f['Curves/Color/Sigmap=0.954'][...]
            curveOffset_f = f['Curves/Color/OffsetFluxDensity'][...]
            curveOffset_mean = f['Curves/Color/OffsetMean'][...]
            curveOffset_1sig = f['Curves/Color/Offsetp=0.683'][...]
            curveOffset_2sig = f['Curves/Color/Offsetp=0.954'][...]
  
            bands = find_bandname(h5file)
            band1 = bands[0] if len(bands) >= 2 else "Unknown"
            band2 = bands[1] if len(bands) >= 2 else "Unknown"

            data = stats2D(band1, band2, modelType, are_ratios, 
                           param_names, knot_positions,
                           sigma_positions, offset_positions, params_init,
                           best_like, params_best, params_mean, params_1sig,
                           curve1D_f, curve1D_mean, curve1D_1sig,
                           curve1D_2sig, curveSigma_f, curveSigma_mean, 
                           curveSigma_1sig, curveSigma_2sig, curveOffset_f, 
                           curveOffset_mean, curveOffset_1sig,
                           curveOffset_2sig)
        else:
            raise ValueError("Unknown model type {0:s}".format(modelType))
    return data


def print_summary(stats):
    """ Prints summary statistics"""

    is2D = stats.ModelType == "numberCountsDoubleLogNormal"

    # Print
    print("#### Model ####")
    print(" Model type: {0:s}".format(stats.ModelType))
    # First -- model knots, which have positions to show
    nknots = len(stats.KnotPositions)
    if is2D:
        nsigmas = len(stats.SigmaKnotPositions)
        noffsets = len(stats.OffsetKnotPositions)
        allknots = np.hstack([stats.KnotPositions, 
                              stats.SigmaKnotPositions,
                              stats.OffsetKnotPositions])
    else:
        allknots = stats.KnotPositions
    ndparams = len(allknots)
    nmult = 2 if is2D else 1
    nparams = len(stats.ParamsBest)

    print("%-11s  %-6s  %5s %6s %6s %6s %6s" %
          ("#Param", "Knot", "Init", "Best", "Fit", "Unc+", "Unc-"))
    for i in range(ndparams):
        print("%-11s  %-6g  %5.2f %6.3f %6.3f %+6.3f %+6.3f" %
              (stats.ParamNames[i], allknots[i], 
               stats.ParamsInit[i], stats.ParamsBest[i],
               stats.ParamsMean[i], stats.Params1Sig[i, 1],
               stats.Params1Sig[i, 0]))

    print("#### Sigma multipliers ####")
    print("%-11s  %5s %6s %6s %6s %6s" %
          ("#Param", "Init", "Best", "Fit", "Unc+", "Unc-"))
    for i in range(ndparams, ndparams + nmult):
        print("%-11s  %5.2f %6.3f %6.3f %+6.3f %+6.3f" %
              (stats.ParamNames[i], stats.ParamsInit[i], 
               stats.ParamsBest[i], stats.ParamsMean[i], 
               stats.Params1Sig[i, 1], stats.Params1Sig[i, 0]))

    print("#### Bonus params ####")
    print("%-11s  %8s %8s %7s %7s" %
          ("#Param", "Best", "Fit", "Unc+", "Unc-"))
    for i in range(ndparams + nmult, nparams):
        print("%-11s  %8.3f %8.3f %+7.3f %+7.3f" %
              (stats.ParamNames[i], stats.ParamsBest[i], 
               stats.ParamsMean[i], stats.Params1Sig[i, 1], 
               stats.Params1Sig[i, 0]))

    print("Best log like: {0:0.2f}".format(stats.BestLikelihood))


def euclideanize(knotpos, knotvals):
    return knotvals + 2.5 * np.log10(knotpos)


def otherdata_helper(other, ax, ms, fmt, label):
    othererr = np.empty((2, other.shape[0]), dtype=np.float)
    othererr[0, :] = np.abs(other[:, 3])  # - errors
    othererr[1, :] = other[:, 2]  # + errors
    ax.errorbar(other[:, 0], other[:, 1], yerr=othererr,
                fmt=fmt, alpha=0.5, ms=ms, mew=1, fillstyle='full',
                label=label)


def band1_plot(ax, stats, errorsnake=None,
               showglenn=False, showbeth=False,
               showoliver=False, showinit=False, euclidean=False,
               skipfirst=False):
    """ Does plot for 1D band model"""
    from pkg_resources import resource_filename

    npointsinterp = 100

    # Read in comparison data if needed
    bandname = stats.Band1 if hasattr(stats, 'Band1') else stats.Band
    if bandname not in ['PSW', 'PMW', 'PLW']:
        if showglenn:
            errstr = "Don't have Gleen et al. 2010 for band %s" % bandname
            raise ValueError(errstr)
        if showoliver:
            errstr = "Don't have Oliver et al. 2010 for band %s" % bandname
            raise ValueError(errstr)
        if showbeth:
            errstr = "Don't have Bethermin et al. 2012 for band %s"
            raise ValueError(errstr % bandname)
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
    sidx = 1 if skipfirst else 0
    n1D = len(stats.KnotPositions)
    knots1Dplot = stats.KnotPositions[sidx:].copy()
    fit1Dplot = stats.ParamsMean[sidx:n1D].copy()
    unc_plus_plot = stats.Params1Sig[sidx:n1D, 1].copy()
    unc_minus_plot = stats.Params1Sig[sidx:n1D, 0].copy()
    initplot = stats.ParamsInit[sidx:n1D].copy()

    errsnake_f = stats.Snake1DFlux.copy()
    errsnake_mn = stats.Snake1DMean.copy()
    errsnake_m1 = stats.Snake1D1Sig[:, 0].copy()
    errsnake_p1 = stats.Snake1D1Sig[:, 1].copy()
    errsnake_m2 = stats.Snake1D2Sig[:, 0].copy()
    errsnake_p2 = stats.Snake1D2Sig[:, 1].copy()

    if euclidean:
        fit1Dplot = euclideanize(knots1Dplot, fit1Dplot)
        initplot = euclideanize(knots1Dplot, initplot)
        if showglenn:
            glenn[:, 1] = euclideanize(glenn[:, 0], glenn[:, 1])
        if showbeth:
            beth[:, 1] = euclideanize(beth[:, 0], beth[:, 1])
        if showoliver:
            oliver[:, 1] = euclideanize(oliver[:, 0], oliver[:, 1])
        errsnake_mn = euclideanize(errsnake_f, errsnake_mn)
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


def sigma_plot(ax, stats):
    n1D = len(stats.KnotPositions)
    nsigma = len(stats.SigmaKnotPositions)
    if nsigma < 2:
        raise ValueError("Don't call sigma_plot with one sigma knot")

    # There is a tricky issue here: ParamsMean are in model
    # space, but if we are plotting ratios then we need them
    # to be in ratio space.  If that is the case, we will just
    # interpolate on the error snake at the knot positions.
    fiterr = np.empty((2, nsigma), dtype=np.float32)
    if stats.AreRatios:
        i1D = interp1d(stats.SnakeSigmaFlux, stats.SnakeSigmaMean,
                       kind='linear')
        sigmavals = i1D(stats.SigmaKnotPositions)
        i1D = interp1d(stats.SnakeSigmaFlux, stats.SnakeSigma1Sig[:, 0])
        fiterr[0, :] = sigmavals - i1D(stats.SigmaKnotPositions)
        i1D = interp1d(stats.SnakeSigmaFlux, stats.SnakeSigma1Sig[:, 1])
        fiterr[1, :] = i1D(stats.SigmaKnotPositions) - sigmavals
    else:
        # Can just take them as they are
        sigmavals = stats.ParamsMean[n1D:(n1D + nsigma)]
        fiterr[0, :] = np.abs(stats.Params1Sig[n1D:(n1D + nsigma), 0])
        fiterr[1, :] = stats.Params1Sig[n1D:(n1D + nsigma), 1]
    ax.plot(stats.SigmaKnotPositions, sigmavals, 'ro', alpha=0.3)
    ax.errorbar(stats.SigmaKnotPositions, sigmavals, 
                yerr=fiterr, fillstyle='full',
                fmt='ro', alpha=0.7, label="This Fit")

    ax.plot(stats.SnakeSigmaFlux, stats.SnakeSigmaMean, 'r', alpha=0.3)
    ax.fill_between(stats.SnakeSigmaFlux, stats.SnakeSigma2Sig[:, 0],
                    stats.SnakeSigma2Sig[:, 1],facecolor='grey', 
                    alpha=0.2, lw=0, edgecolor='grey')
    ax.fill_between(stats.SnakeSigmaFlux, stats.SnakeSigma1Sig[:, 0],
                    stats.SnakeSigma1Sig[:, 1],
                    facecolor='grey', alpha=0.4, lw=0, edgecolor='grey')

    
    ax.set_title("Sigma model")
    ax.set_xlabel("{0:s} Flux Density [Jy]".format(stats.Band1))
    ax.set_xscale('log')
    if stats.AreRatios:
        ax.set_ylabel(r"""$\sigma_{f_2 / f_1}$""")
    else:
        ax.set_ylabel(r"""$\sigma_{\log f_2 / f_1}$""")
    ax.set_xlim(0.5 * stats.SigmaKnotPositions.min(),
                1.5 * stats.SigmaKnotPositions.max())


def offset_plot(ax, stats):
    n1D = len(stats.KnotPositions)
    nsigma = len(stats.SigmaKnotPositions)
    noffset = len(stats.OffsetKnotPositions)
    if noffset < 2:
        raise ValueError("Don't call offset_plot with one offset knot")

    idx1 = n1D + nsigma
    idx2 = idx1 + noffset

    # See discussion in sigma_plot
    fiterr = np.empty((2, noffset), dtype=np.float32)
    if stats.AreRatios:
        i1D = interp1d(stats.SnakeOffsetFlux, stats.SnakeOffsetMean,
                       kind='linear')
        offsetvals = i1D(stats.OffsetKnotPositions)
        i1D = interp1d(stats.SnakeOffsetFlux, stats.SnakeOffset1Sig[:, 0])
        fiterr[0, :] = offsetvals - i1D(stats.OffsetKnotPositions)
        i1D = interp1d(stats.SnakeOffsetFlux, stats.SnakeOffset1Sig[:, 1])
        fiterr[1, :] = i1D(stats.OffsetKnotPositions) - offsetvals
    else:
        # Can just take them as they are
        offsetvals = stats.ParamsMean[idx1:idx2]
        fiterr[0, :] = np.abs(stats.Params1Sig[idx1:idx2, 0])
        fiterr[1, :] = stats.Params1Sig[idx1:idx2, 1]
    ax.plot(stats.OffsetKnotPositions, offsetvals, 'ro', alpha=0.3)
    ax.errorbar(stats.OffsetKnotPositions, offsetvals, 
                yerr=fiterr, fillstyle='full',
                fmt='ro', alpha=0.7, label="This Fit")

    ax.plot(stats.SnakeOffsetFlux, stats.SnakeOffsetMean, 'r', alpha=0.3)
    ax.fill_between(stats.SnakeOffsetFlux, stats.SnakeOffset2Sig[:, 0],
                    stats.SnakeOffset2Sig[:, 1],facecolor='grey', 
                    alpha=0.2, lw=0, edgecolor='grey')
    ax.fill_between(stats.SnakeOffsetFlux, stats.SnakeOffset1Sig[:, 0],
                    stats.SnakeOffset1Sig[:, 1],
                    facecolor='grey', alpha=0.4, lw=0, edgecolor='grey')

    ax.set_title("Offset model")
    ax.set_xlabel("{0:s} Flux Density [Jy]".format(stats.Band1))
    ax.set_xscale('log')
    if stats.AreRatios:
        ax.set_ylabel(r"""$\left<f_2 / f_1\right>$""")
    else:
        ax.set_ylabel(r"""$\left<\log f_2 / f_1\right>$""")
    ax.set_ylabel(r"""$\left<f_2 / f_1\right>$""")
    ax.set_xlim(0.5 * stats.OffsetKnotPositions.min(),
                1.5 * stats.OffsetKnotPositions.max())


def make_plots(stats, showglenn=False, showbeth=False, showoliver=False,
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
    is2D = stats.ModelType == "numberCountsDoubleLogNormal"
    if is2D:
        do_sigmaplot = len(stats.SigmaKnotPositions) > 1
        do_offsetplot = len(stats.OffsetKnotPositions) > 1

    # The band 1 uses 2 panels, others use at most 1
    n_panels = 2 + (1 if do_sigmaplot else 0) + (1 if do_offsetplot else 0)
    figsize = {2: (10, 5), 3: (13, 5), 4: (16, 5)}

    f = plt.figure(figsize=figsize[n_panels])
    gs = gridspec.GridSpec(1, n_panels)
    gs.update(left=0.08, right=0.92, wspace=0.35)

    if is2D:
        # 1D doesn't really need a supertitle
        f.suptitle("{0:s} vs. {1:s}".format(stats.Band1, stats.Band2))

    ax1 = plt.subplot(gs[:, 0:2])
    band1_plot(ax1, stats, showglenn=showglenn, showbeth=showbeth,
               showoliver=showoliver, showinit=showinit,
               euclidean=euclidean, skipfirst=skipfirst)

    if do_sigmaplot:
        ax2 = plt.subplot(gs[:, 2])
        sigma_plot(ax2, stats)

    if do_offsetplot:
        ax3 = plt.subplot(gs[:, 3 if do_sigmaplot else 2])
        offset_plot(ax3, stats)

    return f
