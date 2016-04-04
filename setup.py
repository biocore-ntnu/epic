import sys
from setuptools import setup, find_packages
from Cython.Build import cythonize

install_requires = ["scipy", "pandas", "numpy", "natsort", "ebs", "docopt",
                    "joblib", "pyfaidx"]

if sys.version_info[0] == 2:
    install_requires.append("functools32")

setup(
    name="bioepic",
    packages=find_packages(),
    # ext_modules=cythonize(
    #      "epic/statistics/add_to_island_expectations_cython.pyx"),
    scripts=["bin/epic", "bin/epic-effective"],
    version="0.0.4",
    description="Chip-Seq broad peak/domain finder.",
    author="Endre Bakken Stovner",
    author_email="endrebak@stud.ntnu.no",
    url="http://github.com/endrebak/epic",
    keywords=["ChIP-Seq"],
    license=["MIT"],
    install_requires=install_requires,
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Development Status :: 2 - Pre-Alpha",
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
