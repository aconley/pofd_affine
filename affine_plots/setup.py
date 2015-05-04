from distutils.core import setup

import sys
major, minor1, minor2, release, serial = sys.version_info

if (major < 3) and (minor1 < 7):
    raise SystemExit("affine_plots requires at least python 2.7")

setup(
    name="affine_plots",
    version="0.1.0",
    author="Alexander Conley",
    author_email="alexander.conley@colorado.edu",
    packages=["affine_plots"],
    package_data = {'affine_plots': ['resources/*.txt']},
    scripts = ["affine_plots/make_affine_plots.py"],
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
    requires = ['numpy (>1.7.0)', 'scipy (>0.8.0)', 
                'h5py (>2.0.0)', 'matplotlib (>1.3.0)']
)
