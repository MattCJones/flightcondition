#!/usr/bin/env python
"""
Setup file for flightcondition.
"""

from setuptools import setup, find_packages

# Safely load __version__ without importing package
with open("src/flightcondition/version.py") as f:
    exec(f.read())

if __name__ == "__main__":
    with open('README.rst', 'r') as file:
        long_description_ = file.read()

    print("Packages found:", find_packages(where='src'))  # -v option to see
    setup(
        name='flightcondition',
        version=__version__,
        packages=find_packages(where='src'),
        package_dir={'flightcondition': "src/flightcondition"},
        description=("Airspeed conversions (true/calibrated/equivalent/Mach), "
                     "atmospheric data, and more with built-in unit checking."
                     ),
        long_description=long_description_,
        license="MIT License",
        url="https://github.com/MattCJones/flightcondition",
        author="Matthew C. Jones",
        author_email='matt.c.jones.aoe@gmail.com',
        scripts=None,
        classifiers=[
            "Development Status :: 4 - Beta",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Natural Language :: English",
            "Operating System :: OS Independent",
            "Topic :: Scientific/Engineering",
            "Topic :: Utilities",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            ],
        keywords="utility aerospace engineering design",
        install_requires=('pint', 'numpy'),
    )
