""" 
Assess contacts and determine which sets can be passsed to gibbs.

Thismodule provides the `Assess` class, which determines from the contact file
the residues that have enough data for the Gibbs sampling step. This can also be used 
to identify cases where long contact events occur but there is not sufficient 
sampling to perform Gibbs sampling.
"""
import os
import gc
import pickle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy.random import default_rng
from tqdm import tqdm
from MDAnalysis.analysis.base import Results
from basicrta.util import confidence_interval
import multiprocessing
from multiprocessing import Pool, Lock
import MDAnalysis as mda
from basicrta import istarmap

gc.enable()
mpl.rcParams['pdf.fonttype'] = 42
rng = default_rng()

class Assess(object):
    r"""
    Class used to analyze contact file before sending data to Gibbs samplers.
    This class determines which residues will be passed to the Gibbs sampling
    step, and makes notes of the ones that will not. In particular, it will
    determine possibly important residues for which an insufficient number of 
    unbinding events occurred, such as an observed binding but no observed 
    unbinding.
    """

    def __init__(self, contacts, resid=None, niter=110000):
        self.contact_file = contacts
        self.resid = resid
        self.niter = niter
        
        # load contacts file
        with open(contacts, 'r+b') as f:
            self.contacts = pickle.load(f)

    def collect(self):
        time_sorts = self.contacts[:,-1].argsort()[::-1]
        inds = np.unique(self.contacts[time_sorts,0], return_index=True)[1]

        self.longest_resids = self.contacts[time_sorts, 0][sorted(inds)].astype(int)
        self.longest_times = self.contacts[time_sorts, -1][sorted(inds)]
        self.times = [self.contacts[self.contacts[:,0] == resid, -1] for resid in self.longest_resids]
        self.num_events = np.array([len(time) for time in self.times])

    def get_pose(self, resid, n=0):
        """Get the average pose for a given resid corresponding to the nth 
        longest binding event observed. Note: Longest time is index 0, second 
        longest index 1, etc.
        """
        from MDAnalysis.analysis.rms import rmsd
        from basicrta.util import get_suffix, get_code
 
        top = self.contacts.dtype.metadata['top']
        traj = self.contacts.dtype.metadata['traj']
        dt = self.contacts.dtype.metadata['ts']
        ag1 = self.contacts.dtype.metadata['ag1']
        ag2 = self.contacts.dtype.metadata['ag2']
        cutoff = self.contacts.dtype.metadata['cutoff']

        frames, lipid = self._get_longest_nth_time_info(resid, n)
        ag1_inds = self.contacts.dtype.metadata['ag1'].indices
        ag2_sel_inds = ag2.select_atoms(f'resid {lipid}').indices
        
        u = mda.Universe(f'{top}', f'{traj}')
        traj_len = len(frames)
    
        # redefine selections
        new_ag1 = u.atoms[ag1_inds]
        new_ag2 = u.atoms[ag2_sel_inds]

        # get average position
        positions = 0
        for i, ts in enumerate(u.trajectory[frames]):
            positions += new_ag2.positions
        positions = positions / traj_len
        
        # get closest frame to average
        rmsds = np.empty(traj_len)
        for i, ts in enumerate(u.trajectory[frames]):
            rmsds[i] = rmsd(new_ag2.positions, positions)
        
        min_index = rmsds.argmin()
        write_frame = frames[min_index]
        u.trajectory[write_frame]
        
        reslet = get_code(ag1.select_atoms(f'resid {resid}').residue.resname)
        residue = f'{reslet}{resid}'
        #outdir = f'basicrta-{cutoff}/{residue}/'
        out_name = f'{residue}_{n+1}{get_suffix(n+1)}_longest.pdb'
        (new_ag1 +new_ag2).atoms.write(out_name)

    def _get_longest_nth_time_info(self, resid, n):
        dt = self.contacts.dtype.metadata['ts']
        tmp = self.contacts[self.contacts[:,0] == resid] 
        sorted_tmp = tmp[tmp[:,-1].argsort()[::-1]]
        start_frame = round(sorted_tmp[n, 2] / dt)
        end_frame = round(start_frame + sorted_tmp[n, 3] / dt)
        lipid = int(sorted_tmp[n, 1])
        return np.array(range(start_frame, end_frame)), lipid 
    
    #def compare_poses(self):
    #    for resid in self.longest_resids:
    #        if not os.path.exists()
    #        self.get_pose(resid)

def get_parser():
    import argparse
    parser = argparse.ArgumentParser(description="""run gibbs samplers for all
                                     or a specified residue present in the
                                     contact map""",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group('required arguments')

    required.add_argument('--contacts', required=True, help="""Contact file
                          produced from `basicrta contacts`, default is
                          contacts_{cutoff}.pkl""")
    parser.add_argument('--resid', type=int, help="""run gibbs sampler for
                          this residue. Will collect cutoff from contact file
                          name.""")
    parser.add_argument('--niter', type=int, default=110000, help="""number of
                          iterations to use for the gibbs sampler""")
    # this is to make the cli work, should be just a temporary solution
    parser.add_argument('assses', nargs='?', help=argparse.SUPPRESS)
    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()

    contact_path = os.path.abspath(args.contacts)
    cutoff = args.contacts.split('/')[-1].strip('.pkl').split('_')[-1]

    ParallelGibbs(contact_path, nproc=args.nproc, ncomp=args.ncomp,
                  niter=args.niter).run(run_resids=args.resid)

if __name__ == '__main__':
    exit(main())
