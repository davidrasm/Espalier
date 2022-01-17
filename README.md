# Espalier

Espalier is a Python package for working with discordant phylogenetic trees using maximum agreement forests. The package can be used to compute maximum agreement forests, compute subtree-prune-and-regraft (SPR) distances, reconcile discordant trees and reconstruct approximate ancestral recombination graphs (ARGs).

For further information on how to use Espalier, please see the [primer][primer] and [documentation][docs] pages.

[primer]: <https://espalier.readthedocs.io/en/latest/primer.html>
[docs]: <https://espalier.readthedocs.io/en/latest/>

For more information on the algorithms behind Espalier please see:

*Rasmussen, D.A. and Guo, F. Espalier: Efficient tree reconciliation and ARG reconstruction using maximum agreement forests. 2022.*

## Installation and set up

The easiest way to install Espalier is through [pip][pip]:
```
$ pip install Espalier
```

pip will install the other Python packages required by Espalier. 

However, Espalier also requires but is not packaged with RAxML-NG. Instructions for installing RAxML-NG can be found [here][raxml-ng].

If Espalier was installed correctly, the command
```
$ espalier --help
```

should return the options for Espalier's command line interface.

[pip]: <https://pypi.org/project/pip/>
[raxml-ng]: <https://github.com/amkozlov/raxml-ng>