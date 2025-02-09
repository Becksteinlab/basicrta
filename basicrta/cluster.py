import os
import gc
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool, Lock
import MDAnalysis as mda
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
    """
    
    def __init__(self, niter, prot, cutoff, gskip):
        self.residues = {}
        self.niter = niter
        self.prot = prot
        self.cutoff = cutoff
        self.gskip = gskip

    def __getitem__(self, item):
        return getattr(self, item)

    def _single_residue(self, adir, process=False):
        if os.path.exists(f'{adir}/gibbs_{self.niter}.pkl'):
            try:
                result = f'{adir}/gibbs_{self.niter}.pkl'
                g = Gibbs().load(result)
                if process:
                    g.gskip = self.gskip
                    g.process_gibbs()
            except ValueError:
                result = None
        else:
            print(f'results for {adir} do not exist')
            result = None
        return result

    def _single_result(self, adir):
        result = self._single_residue(adir)
        residue = adir.split('/')[-1]
        self.residues[residue] = result


    def reprocess(self, nproc=1):
        """Rerun processing and clustering on :class:`Gibbs` data.

        :param nproc: Number of processes to use in clustering results for all
                      residues.
        :type nproc: int
        """
        from glob import glob

        dirs = np.array(glob(f'basicrta-{self.cutoff}/?[0-9]*'))
        sorted_inds = (np.array([int(adir.split('/')[-1][1:]) for adir in dirs])
                       .argsort())
        dirs = dirs[sorted_inds]
        inarr = np.array([[adir, True] for adir in dirs])
        with (Pool(nproc, initializer=tqdm.set_lock,
                   initargs=(Lock(),)) as p):
            try:
                for _ in tqdm(p.istarmap(self._single_residue, inarr),
                              total=len(dirs), position=0,
                              desc='overall progress'):
                    pass
            except KeyboardInterrupt:
                pass

    def collect_results(self, nproc=1):
        """Collect names of results for each residue in the `basicrta-{cutoff}`
        directory in a dictionary stored in :attr:`ProcessProtein.results`.
        """
        from glob import glob

        dirs = np.array(glob(f'basicrta-{self.cutoff}/?[0-9]*'))
        sorted_inds = (np.array([int(adir.split('/')[-1][1:]) for adir in dirs])
                       .argsort())
        dirs = dirs[sorted_inds]
        with (Pool(nproc, initializer=tqdm.set_lock,
                   initargs=(Lock(),)) as p):
            try:
                for _ in tqdm(p.istarmap(self._single_result, dirs),
                              total=len(dirs), position=0,
                              desc='overall progress'):
                    pass
                #for adir in tqdm(dirs, desc='collecting results'):
                #    result = self._single_residue(adir)
                #    residue = adir.split('/')[-1]
                #    self.residues[residue] = result
                except KeyboardInterrupt:
                    pass

    def get_taus(self):
        r"""Get :math:`\tau` and 95\% confidence interval bounds for the slowest
        process for each residue. 
        
        :returns: Returns a tuple of the form (:math:`\tau`, [CI lower bound, 
                  CI upper bound])
        :rtype: tuple
        
        """
        from basicrta.util import get_bars

        taus = []
        for res in tqdm(self.residues, total=len(self.residues)):
            if self.residues[res] is None:
                result = [0, 0, 0]
            else:
                try:
                    gib = Gibbs().load(self.residues[res])
                    result = gib.estimate_tau()
                except AttributeError:
                    result = [0, 0, 0]
            taus.append(result)
        taus = np.array(taus)
        bars = get_bars(taus)
        return taus[:, 1], bars

    def write_data(self, fname='tausout'):
        r"""Write :math:`\tau` values with 95\% confidence interval to a numpy
        file with the format [`sel1` resid, :math:`\tau`, CI lower bound, CI
        upper bound].

        :param fname: Filename to save data to.
        :type fname: str, optional
        """
        taus, bars = self.get_taus()
        keys = self.residues.keys()
        residues = np.array([int(res[1:]) for res in keys])
        data = np.stack((residues, taus, bars[0], bars[1]))
        np.save(fname, data.T)

    def plot_protein(self, **kwargs):
        r"""Plot :math:`\tau` vs resid. kwargs are passed to the 
        :meth:`plot_protein` method of `util.py`. These can be used to change
        the labeling cutoff, y-limit of the plot, scale the figure, and set
        major and minor ticks.
        """
        from basicrta.util import plot_protein
        if len(self.residues) == 0:
            print('run `collect_residues` then rerun')

        taus, bars = self.get_taus()
        residues = list(self.residues.keys())
        residues = [res.split('/')[-1] for res in residues]

        exclude_inds = np.where(bars < 0)[1] 
        taus = np.delete(taus, exclude_inds)
        bars = np.delete(bars, exclude_inds, axis=1)
        residues = np.delete(residues, exclude_inds)

        plot_protein(residues, taus, bars, self.prot, **kwargs)

    def b_color_structure(self, structure):
        r"""Add :math:`\tau` to b-factors in the specified structure. Saves
        structure with b-factors to `tau_bcolored.pdb`. 
        """
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


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--nproc', type=int, default=1)
    parser.add_argument('--cutoff', type=float)
    parser.add_argument('--niter', type=int, default=110000)
    parser.add_argument('--prot', type=str, default=None, nargs='?')
    parser.add_argument('--label-cutoff', type=float, default=3,
            dest='label_cutoff',
            help='Only label residues with tau > '
            'LABEL-CUTOFF * <tau>. ')

    parser.add_argument('--structure', type=str, nargs='?')
    args = parser.parse_args()

    pp = ProcessProtein(args.niter, args.prot, args.cutoff)
    pp.reprocess(nproc=args.nproc)
    pp.collect_results()
    pp.write_data()
    pp.plot_protein(label_cutoff=args.label_cutoff)
