.. basicrta documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to basicrta's documentation!
=========================================================
``basicrtakit`` is a package to analyze binding kinetics in molecular dynamics 
simulations. The analysis uses an exponential mixture model and Bayesian 
nonparametric inference to determine time-scales of the underlying binding 
processes. The output is processed to give frames belonging to a given process 
with an associated time-scale, which can be further analyzed.

This package was written to perform the analyses in [Sexton_2024]_. 


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   getting_started
   tutorial
   api


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. [Sexton_2024] Sexton, R.; Fazel, M.; Schweiger, M.; Pressé, S.; 
   Beckstein, O. Bayesian Nonparametric Analysis of Residence Times for 
   Protein-Lipid Interactions in Molecular Dynamics Simulations. bioRxiv. November 9, 
   2024. https://doi.org/10.1101/2024.11.07.622502.

