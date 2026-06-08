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
URL = 'https://github.com/davidrasm/Espalier'

LICENSE = 'Apache License 2.0'
DESCRIPTION = 'Espalier is a Python package for working with discordant phylogenetic trees using maximum agreement forests'
LONG_DESCRIPTION = (HERE / "README.md").read_text()
LONG_DESC_TYPE = "text/markdown"

INSTALL_REQUIRES = [
      'numpy>=2.4.6',
      'DendroPy>=5.0.8',
      'pandas>=3.0.3',
      'tskit>=1.0.3',
      'biopython>=1.87',
      'scipy>=1.17.1',
      'click>=8.4.1',
]

EXTRAS_REQUIRE = {
      'simulation': [
            'msprime>=1.4.2',
            'pyvolve>=1.1.0',
      ],
      'dev': [
            'msprime>=1.4.2',
            'pyvolve>=1.1.0',
            'pytest>=9.0.3',
      ],
}

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
      extras_require=EXTRAS_REQUIRE,
      python_requires='>=3.11',
      #packages=['Espalier'],
      packages=find_packages(),
      entry_points={
            'console_scripts': [
                  'espalier = Espalier.espalier:main',
            ],
            }
      )
