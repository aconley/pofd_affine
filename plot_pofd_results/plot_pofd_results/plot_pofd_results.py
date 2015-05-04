#!python

from __future__ import print_function

from plot_affine_results import *

if __name__ == "__main__":
    import argparse
    import plot_affine_results

    desc = """Plot and print the results of a P(D) fit"""

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
        outfile = os.path.splitext(parse_results.h5file)[0] + '.pdf'
        f.savefig(outfile)
    else:
        f.savefig(outfile)
