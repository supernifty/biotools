#!/usr/bin/env python

import setuptools
import pathlib

name = "biotools"
version = "0.1"
release = "0.1.0"
here = pathlib.Path(__file__).parent.resolve()

setuptools.setup(
    name=name,
    version=version,
    packages=setuptools.find_packages()
)
