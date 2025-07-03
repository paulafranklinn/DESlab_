#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Setup script for DESlab
"""

from setuptools import setup, find_packages

deslab_email = "deslab@ufrj.br"

deslab_description = "DESlab is a scientific computing program for discrete event systems (DES) modeled as automata."

deslab_long_description = """
DESlab
========
DESlab is a scientific computing program written in Python
for the development of algorithms for analysis and synthesis
of discrete event systems (DES) modeled as automata.
It integrates automata, graph algorithms, and numerical
calculations. DESlab also allows the definition of symbolic
variables and incorporates concise instructions to manipulate,
analyze and visualize them, with syntax close to DES theory.
"""

setup(
    name="DESlab",
    version="0.0.4",
    maintainer="Leonardo Clavijo, Daniel Ramos Garcia, Joao Carlos Basilio, Lilian Kawakami",
    maintainer_email=deslab_email,
    author="Leonardo Clavijo, Daniel Ramos Garcia, Joao Carlos Basilio, Lilian Kawakami",
    author_email=deslab_email,
    description=deslab_description,
    long_description=deslab_long_description,
    long_description_content_type="text/plain",
    license="BSD",
    platforms=['any'],
    url="http://www.dee.ufrj.br/lca/",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'deslab.graphics': ['output/empty.dat', 'working/*.py'],
        'deslab': ['docs/*.pdf'],
    },
    install_requires=[
        "networkx>=2.8",
        "pandas>=2.0",
        "pyparsing>=3.0",
        "pydot>=1.4",
        "portion>=2.0"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)
