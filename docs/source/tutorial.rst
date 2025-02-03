.. |AA| unicode:: U+212B 

========
Tutorial
========

Contact Analysis
================ 

The basicrta workflow starts with collecting contacts between two atom groups
``sel1`` and ``sel2`` based on a single ``cutoff`` using contacts.py.::
  python -m basicrta.contacts --top top.pdb --traj traj.xtc --sel1 "protein"
  --sel2 "resname CHOL" --cutoff 7.0 --nproc 4

This will create two contact maps: ``contacts.npy`` and ``contacts_7.0.pkl``.
The first is used to map all contacts defined with a 10 |AA| cutoff, upon which
the cutoff is imposed and a new contact map is created for each desired cutoff
value (for storage reasons this may change in the future). The protein residues
in the topology should be correctly numbered, as these will be the names of
directories where results are stored.  

Gibbs Sampler
=============

Next the contact map is used to collect the contacts between a specified residue
of ``sel1`` and each residue of the ``sel2`` group. The contact durations
(residence times) are then used as input data for the Gibbs sampler. A specific
residue of ``sel1`` can be used to run a Gibbs sampler for only that residue :: 
  python -m basicrta.gibbs --contacts contacts_7.0.pkl --nproc 5 --resid 313
or if ``resid`` is left out, a Gibbs sampler will be executed for all ``sel1``
residues in the contact map. ::
  python -m basicrta.gibbs --contacts contacts_7.0.pkl --nproc 5


Clustering
==========

Next the samples obtained from the Gibbs sampler are processed and clustered. 
::
  python -m basicrta.cluster --niter 110000 --nproc 3 --cutoff 7.0 --prot b2ar

The ``prot`` argument is used to create rectangles in the :math:`\tau` vs resid
plot that correspond to the TM segments of the protein (figures 7-10). Your
protein can be added to ``basicrta/basicrta/data/tm_dict.txt`` in the same
format as the existing proteins. 

``basicrta.cluster`` will process the Gibbs samplers, compute :math:`\tau` for
each residue, plot :math:`\tau` vs resid, and write the data to ``tausout.npy``
(see next section). If a structure is passed to the script, the b-factors of the
residues will be populated with the appropriate :math:`\tau`.

Rate/Tau Estimates
==================

Estimates for :math:`\tau` are obtained by using the ``write_data()`` method of
the ``ProcessProtein`` class contained in ``cluster.py``. Data is saved as a ``.npy``
file and contains [protein residue, tau, CI lower bound, CI upper bound] for
each residue of ``sel1`` analyzed (given dataset size is sufficiently large
(~50 points)). Writing out all model parameters to a ``.npy`` file is currently
not implemented, but will be possible in later versions.


Kinetic Mapping
===============

The kinetically mapped trajectory and weighted densities can be created using 
``kinetics.py``. ::
  python -m basicrta.kinetics --gibbs basicrta_7.0/W313/gibbs_110000.pkl
  --contacts contacts_7.0.pkl --wdensity

To create only the mapped trajectory, leave out the ``wdensity`` flag.  
::
  python -m basicrta.kinetics --gibbs basicrta_7.0/W313/gibbs_110000.pkl
  --contacts contacts_7.0.pkl

This can also be done using the ``top_n`` most likely frames belonging to each
component of the exponential mixture model. ::
  python -m basicrta.kinetics --gibbs basicrta_7.0/W313/gibbs_110000.pkl
  --contacts contacts_7.0.pkl --top_n 500

Weighted densities can be computed over the ``top_n`` frames, over the whole
trajectory, or by using the ``step`` argument in combination with ``top_n`` or
the whole trajectory. ::
  python -m basicrta.kinetics --gibbs basicrta_7.0/W313/gibbs_110000.pkl
  --contacts contacts_7.0.pkl --step 100

Supplemental Scripts
====================

Slurm scripts for submitting Gibbs sampler jobs to distributed systems are
located in the ``basicrta/scripts`` directory. Note that some of these were
rewritten and may contain slight errors, testing still needs to be done (but
can be used with some slight modifications).
