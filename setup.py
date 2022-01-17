"""
Created on Fri Nov 26 11:07:53 2021

@author: david
"""
import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent

VERSION = '0.1.1'
PACKAGE_NAME = 'Espalier'
AUTHOR = 'David Rasmussen'
AUTHOR_EMAIL = 'drasmus@ncsu.edu'
URL = 'https://github.ncsu.edu/phylodynamics-lab/Espalier'

LICENSE = 'Apache License 2.0'
DESCRIPTION = 'Espalier reconciles discordant trees and reconcstructs ARGs using maximum agreement forests.'
LONG_DESCRIPTION = (HERE / "README.md").read_text()
LONG_DESC_TYPE = "text/markdown"

# To add for testing???: msprime, pyvolve
INSTALL_REQUIRES = [
      'numpy',
      'dendropy',
      'pandas',
      'tskit',
      'biopython',
      'scipy',
      'click'
]

setup(name=PACKAGE_NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      long_description_content_type=LONG_DESC_TYPE,
      author=AUTHOR,
      license=LICENSE,
      author_email=AUTHOR_EMAIL,
      url=URL,
      install_requires=INSTALL_REQUIRES,
      #packages=['Espalier'],
      packages=find_packages(),
      entry_points={
            'console_scripts': [
                  'espalier = Espalier.espalier:main',
            ],
            }
      )