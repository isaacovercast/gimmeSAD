#!/usr/bin/env python2.7

from setuptools import setup, find_packages
import glob
import re

def requires():
    """ gets packages from requirements.txt """
    with open('requirements.txt') as infile:
        return infile.read().splitlines()

CUR_VERSION = 1.0
setup(
    name="gimmeSAD",
    version=CUR_VERSION,
    url="https://github.com/isaacovercast/gimmeSAD",
    author="Isaac Overcast",
    author_email="iovercast@gc.cuny.edu",
    description="Joint modelling of abundance and genetic diversity. ",
    long_description=open('README.md').read(),
    packages=find_packages(),    
    install_requires=requires(),
    #dependencies=dependency_links(),
    entry_points={
            'console_scripts': [
                'gimmeSAD = gimmeSAD.gimmeSAD:main',
            ],
    },
    license='CC-BY-SA',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
    ],
)
