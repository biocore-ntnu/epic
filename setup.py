import os
import sys
from setuptools import setup, find_packages
# from Cython.Build import cythonize

from epic.version import __version__

install_requires = ["scipy", "pandas>=0.23.0", "numpy", "natsort", "joblib", "pyfaidx", "typing"]

try:
    os.getenv("TRAVIS")
    install_requires.append("coveralls")
except:
    pass

if sys.version_info[0] == 2:
    install_requires.append("functools32")

setup(
    name="bioepic",
    packages=find_packages(),

    scripts=["bin/epic", "bin/epic-effective", "bin/epic-overlaps", "bin/epic-merge", "bin/epic-cluster", "bin/epic-count", "bin/epic-blacklist"],
    package_data={'epic': ['scripts/effective_sizes/*.txt',
                           'scripts/chromsizes/*chromsizes',
                           'scripts/genome.snakefile']},
    version=__version__,
    description="Chip-Seq broad peak/domain finder.",
    author="Endre Bakken Stovner",
    author_email="endrebak85@gmail.com",
    url="http://github.com/endrebak/epic",
    keywords=["ChIP-Seq"],
    license=["MIT"],
    install_requires=install_requires,
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta",
        "Environment :: Other Environment", "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Topic :: Scientific/Engineering"
    ],
    long_description=
    ("Chip-Seq broad peak/domain finder based on SICER. See the url for more info."
     ))
