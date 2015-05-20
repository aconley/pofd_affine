#!python

from __future__ import print_function

if __name__ == "__main__":
    import argparse
    import os.path
    from affine_plots import *
    
    desc = """Plot and print the results of a P(D) fit"""

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('h5file', action='store',
                        help="HDF5 stats file to load")
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

    stats = read_data(parse_results.h5file)

    print_summary(stats)

    if not parse_results.noplot:
      f = make_plots(stats, showglenn=parse_results.glenn,
                     showbeth=parse_results.beth,
                     showoliver=parse_results.oliver,
                     showinit=parse_results.plotinit,
                     euclidean=parse_results.euclidean,
                     skipfirst=parse_results.skipfirst)

      if parse_results.outfile is None:
          outfile = os.path.splitext(parse_results.h5file)[0] + '.pdf'
      else:
          outfile = parse_results.outfile

      print("Writing plots to {0:s}".format(outfile))
      f.savefig(outfile)
