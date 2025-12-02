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
        self.contacts = contacts
        self.resid = resid
        self.niter = niter

    def collect(self):
        # load contacts file
        with open(self.contacts, 'r+b') as f:
            self._contacts = pickle.load(f)

        time_sorts = contacts[:,-1].argsort()[::-1]
        inds = np.unique(contacts[time_sorts,0], return_index=True)[1]

        self.longest_resids = contacts[time_sorts, 0][sorted(inds)].astype(int)
        self.longest_times = contacts[time_sorts, -1][sorted(inds)]
        self.times = [contacts[contacts[:,0] == resid, -1] for resid in self.longest_resids]
        self.num_events = np.array([len(time) for time in self.times])

    def get_pose(self, resid):
        u = MDA.Universe()
    

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
