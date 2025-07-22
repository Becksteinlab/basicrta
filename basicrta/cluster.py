import os
import gc
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool, Lock
from glob import glob
import MDAnalysis as mda
from MDAnalysis.analysis.base import Results
from basicrta import istarmap
from basicrta.gibbs import Gibbs
gc.enable()

"""This module provides the ProcessProtein class, which collects and processes
Gibbs sampler data.
"""

class ProcessProtein(object):
    r"""ProcessProtein is the class that collects and processes Gibbs sampler
    data. This class collects results for all residues in the 
    `basicrta-{cutoff}` directory and can write out a :math:`\tau` vs resid
    numpy array or plot :math:`\tau` vs resid. If a structure is provided,
    :math:`\tau` will be written as b-factors for visualization.

    :param niter: Number of iterations used in the Gibbs samplers
    :type niter: int
    :param prot: Name of protein in `tm_dict.txt`, used to draw TM bars in
                 :math:`tau` vs resid plot.
    :type prot: str, optional
    :param cutoff: Cutoff used in contact analysis.
    :type cutoff: float
    :param gskip: Gibbs skip parameter for decorrelated samples;
                  default from https://pubs.acs.org/doi/10.1021/acs.jctc.4c01522
    :type gskip: int
    :param burnin: Burn-in parameter, drop first N samples as equilibration;
                   default from https://pubs.acs.org/doi/10.1021/acs.jctc.4c01522
    :type burnin: int
    """
    
    def __init__(self, niter, prot, cutoff, 
                 gskip=1000, burnin=10000, 
                 taus=None, bars=None):
        self.residues = Results() # TODO: double-check that we need to use this, it gets set in reprocess/get_taus
        self.niter = niter
        self.prot = prot
        self.cutoff = cutoff
        self.gskip = gskip
        self.burnin = burnin
        self.taus = taus
        self.bars = bars

    def __getitem__(self, item):
        return getattr(self, item)

    def _single_residue(self, adir, process=False):
        if os.path.exists(f'{adir}/gibbs_{self.niter}.pkl'):
            result = f'{adir}/gibbs_{self.niter}.pkl'
            try:
                result = f'{adir}/gibbs_{self.niter}.pkl'
                g = Gibbs().load(result)
                if process:
                    g.gskip = self.gskip
                    g.burnin = self.burnin
                    g.process_gibbs()
                tau = g.estimate_tau()
            except:
                result = None
                tau = [0, 0, 0]
        else:
            result = None
            tau = [0, 0, 0]

        residue = adir.split('/')[-1]
        return residue, tau, result
        #setattr(self.residues, f'{residue}', Results())
        #setattr(self.residues[f'{residue}'], 'file', result)
        #setattr(self.residues[f'{residue}'], 'tau', tau)
        #return self

    def reprocess(self, nproc=1):
        """Rerun processing and clustering on :class:`Gibbs` data.

        :param nproc: Number of processes to use in clustering results for all
                      residues.
        :type nproc: int
        """
        from basicrta.util import get_bars

        dirs = np.array(glob(f'basicrta-{self.cutoff}/?[0-9]*'))
        sorted_inds = (np.array([int(adir.split('/')[-1][1:]) for adir in dirs])
                       .argsort())
        dirs = dirs[sorted_inds]
        inarr = np.array([[adir, True] for adir in dirs])
        with (Pool(nproc, initializer=tqdm.set_lock,
                   initargs=(Lock(),)) as p):
            try:
                residues, taus, results = [], [], []
                for residue, tau, result in tqdm(p.istarmap(self._single_residue, inarr),
                              total=len(dirs), position=0,
                              desc='overall progress'):
                    residues.append(residue)
                    taus.append(tau)
                    results.append(result)
                    gc.collect()
                    pass
            except KeyboardInterrupt:
                pass

        taus = np.array(taus)
        bars = get_bars(taus)
        setattr(self, 'taus', taus[:, 1])
        setattr(self, 'bars', bars)
        setattr(self, 'residues', np.array(residues))
        setattr(self, 'files', np.array(results))

    def get_taus(self, nproc=1):
        r"""Get :math:`\tau` and 95\% confidence interval bounds for the slowest
        process for each residue. 
        
        :returns: Returns a tuple of the form (:math:`\tau`, [CI lower bound, 
                  CI upper bound])
        :rtype: tuple
        
        """
        from basicrta.util import get_bars

        dirs = np.array(glob(f'basicrta-{self.cutoff}/?[0-9]*'))
        sorted_inds = (np.array([int(adir.split('/')[-1][1:]) for adir in dirs])
                       .argsort())
        dirs = dirs[sorted_inds]
        with (Pool(nproc, initializer=tqdm.set_lock,
                   initargs=(Lock(),)) as p):
            try:
                residues, taus, results = [], [], []
                for residue, tau, result in tqdm(p.imap(self._single_residue, dirs),
                                                 total=len(dirs), position=0,
                                                 desc='overall progress'):
                    residues.append(residue)
                    taus.append(tau)
                    results.append(result)
            except KeyboardInterrupt:
                pass
        
        #taus = []
        #for res in tqdm(self.residues, total=len(self.residues)):
        #    taus.append(res.tau)
        
        taus = np.array(taus)
        bars = get_bars(taus)
        setattr(self, 'taus', taus[:, 1])
        setattr(self, 'bars', bars)
        setattr(self, 'residues', np.array(residues))
        setattr(self, 'files', np.array(results))
        return taus[:, 1], bars

    def write_data(self, fname='tausout'):
        r"""Write :math:`\tau` values with 95\% confidence interval to a numpy
        file with the format [`sel1` resid, :math:`\tau`, CI lower bound, CI
        upper bound].

        :param fname: Filename to save data to.
        :type fname: str, optional
        """
        if self.taus is None:
            taus, bars = self.get_taus()

        # Handle residues as numpy array (from reprocess/get_taus methods)
        # TODO: double-check that we need to use res[1:] and can't get this easier
        residues = np.array([int(res[1:]) for res in self.residues])
        data = np.stack((residues, self.taus, self.bars[0], self.bars[1]))
        np.save(fname, data.T)

    def plot_protein(self, **kwargs):
        r"""Plot :math:`\tau` vs resid. kwargs are passed to the 
        :meth:`plot_protein` method of `util.py`. These can be used to change
        the labeling cutoff, y-limit of the plot, scale the figure, and set
        major and minor ticks.
        """
        from basicrta.util import plot_protein

        if self.taus is None:
            self.get_taus()

        residues = self.residues
        residues = [res.split('/')[-1] for res in residues]

        exclude_inds = np.where(self.bars < 0)[1] 
        taus = np.delete(self.taus, exclude_inds)
        bars = np.delete(self.bars, exclude_inds, axis=1)
        residues = np.delete(residues, exclude_inds)

        plot_protein(residues, taus, bars, self.prot, **kwargs)

    def b_color_structure(self, structure):
        r"""Add :math:`\tau` to b-factors in the specified structure. Saves
        structure with b-factors to `tau_bcolored.pdb`. 
        """
        if self.taus is None:
            taus, bars = self.get_taus()

        cis = bars[1]+bars[0]
        errs = taus/cis
        errs[errs != errs] = 0
        residues = list(self.residues.keys())
        u = mda.Universe(structure)

        u.add_TopologyAttr('tempfactors')
        u.add_TopologyAttr('occupancies')
        for tau, err, residue in tqdm(zip(taus, errs, residues)):
            res = u.select_atoms(f'protein and resid {residue[1:]}')
            res.tempfactors = np.round(tau, 2)
            res.occupancies = np.round(err, 2)

        u.select_atoms('protein').write('tau_bcolored.pdb')


if __name__ == "__main__":  #pragma: no cover
    # the script is tested in the test_cluster.py but cannot be accounted for
    # in the coverage report
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--nproc', type=int, default=1)
    parser.add_argument('--cutoff', type=float)
    parser.add_argument('--niter', type=int, default=110000)
    parser.add_argument('--prot', type=str, default=None, nargs='?')
    parser.add_argument('--label-cutoff', type=float, default=3,
            dest='label_cutoff',
            help='Only label residues with tau > '
            'LABEL-CUTOFF * <tau>. ')
    parser.add_argument('--structure', type=str, nargs='?')
    # use  for default values
    parser.add_argument('--gskip', type=int, default=1000, 
                        help='Gibbs skip parameter for decorrelated samples;'
                        'default from https://pubs.acs.org/doi/10.1021/acs.jctc.4c01522')
    parser.add_argument('--burnin', type=int, default=10000, 
                        help='Burn-in parameter, drop first N samples as equilibration;'
                        'default from https://pubs.acs.org/doi/10.1021/acs.jctc.4c01522')

    args = parser.parse_args()

    pp = ProcessProtein(args.niter, args.prot, args.cutoff, 
                        gskip=args.gskip, burnin=args.burnin)
    pp.reprocess(nproc=args.nproc)
    pp.get_taus()
    pp.write_data()
    pp.plot_protein(label_cutoff=args.label_cutoff)
