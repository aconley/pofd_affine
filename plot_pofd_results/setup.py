from distutils.core import setup

import sys
major, minor1, minor2, release, serial = sys.version_info

if (major < 3) and (minor1 < 7):
    raise SystemExit("analyze_coverage requires at least python 2.7")

setup(
    name="plot_pofd_results",
    version="0.1.0",
    author="Alexander Conley",
    author_email="alexander.conley@colorado.edu",
    packages=["plot_pofd_results"],
    package_data = {'plot_pofd_results': ['resources/*.txt']},
    scripts=["plot_pofd_results/plot_pofd_results.py"],
    license="GPL",
    description="Plot results from pofd_affine",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
    requires = ['numpy (>1.7.0)', 'scipy (>0.8.0)', 'h5py (>2.0.0)']
)

