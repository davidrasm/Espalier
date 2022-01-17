##########################
##
## Espalier: A Python package for tree reconciliation and reconstructing ARGs using maximum agreement forests.
##
## Copyright 2021-2022 David A. Rasmussen (drasmus@ncsu.edu)
##
## If using Espalier or this code, please cite:
##
##      Rasmussen, D.A. and Guo, F. Espalier: Efficient tree reconciliation and ARG reconstruction using maximum agreement forests. 2022.
##
############################

import logging

# Setup logging - not sure if __init__ is the best place for this
loglevel = 'INFO'
numeric_level = getattr(logging, loglevel.upper())
logging.basicConfig(format='%(levelname)s: %(message)s', level=numeric_level)
