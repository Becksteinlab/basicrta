from tqdm import tqdm
from basicrta.util import get_start_stop_frames
import numpy as np
import MDAnalysis as mda
import os
import pickle

class CreateTraj(object):
    def __init__(self, contacts):
        self.contacts = contacts

    def create_traj(self, skip=1):
        from numpy.lib.format import open_memmap
        with open(self.contacts, 'rb') as f:
            contacts = pickle.load(f)

        top = contacts.dtype.metadata['top']
        traj = contacts.dtype.metadata['traj']
        ag1 = contacts.dtype.metadata['ag1']
        ag2 = contacts.dtype.metadata['ag2']
        dt = contacts.dtype.metadata['ts']

        nframes = np.round(contacts[:,3]/dt).astype(int)
        start_frames = np.round(contacts[:,2]/dt).astype(int)
        lipinds = contacts[:,1].astype(int)
        del contacts

        u = mda.Universe(top, traj)
        u_reduced = ag1 + ag2.residues[0].atoms
        if not os.path.exists('reduced.pdb'):
            u_reduced.atoms.write('reduced.pdb')

        if not os.path.exists('data.npy'):
            tmplens = [len(np.arange(start, start + N)) for start, N in
                       zip(start_frames, nframes)]
            totlen = sum(tmplens)
            data = open_memmap('data.npy', mode='w+', dtype=np.int64,
                               shape=(totlen, 2))
            q = 0
            for start, N, lipid in tqdm(zip(start_frames, nframes, lipinds),
                                        total=len(nframes), 
                                        desc='writing data'):
                frames = list(range(start, start + N))
                lips = [lipid]*len(frames)
                data[q:q+len(frames), 0] = frames 
                data[q:q+len(frames), 1] = lips
                q += len(frames)
        else:
            data = np.load('data.npy', mmap_mode='r')

        #frames = np.array(frames).astype(int)
        #lips = np.array(lips).astype(int)

        with mda.Writer('full_traj.xtc', len(u_reduced)) as w:
            for i,ts in tqdm(enumerate(u.trajectory[data[:,0]][::skip]),
                             desc='writing trajectory', total=len(data[::skip])):
                w.write((ag1 + ag2.select_atoms(f'resid {data[i*skip, 1]}')).atoms)

#        with mda.Writer('full_traj.xtc', len(u_reduced)) as w:
#            for start, N, lipid in tqdm(zip(start_frames, nframes, lipinds),
#                                        total=len(nframes)):
#                frames = np.arange(start, start + N)
#                for ts in u.trajectory[frames][::skip]:
#                    w.write((ag1 + ag2.select_atoms(f'resid {lipid}')).atoms)

        #frames = [np.arange(start, start + N) for start, N in 
        #          zip(start_frames, nframes)]

        #bframes, eframes = get_start_stop_frames(trajtimes, times, self.ts)
        #tmplens = [len(np.arange(b, e)) for b, e in zip(bframes, eframes)]
        #totlen = sum(tmplens)
        #write_data = open_memmap(self.dataname, mode='w+', dtype=np.float64,
        #                         shape=(totlen, ncomp+2))

        #j = 0
        #for b, e, l, i in tqdm(zip(bframes, eframes, lipinds, indicators),
        #                       total=len(bframes)):
        #    tmp = np.arange(b, e)
        #    tmpl = np.ones_like(np.arange(b, e)) * l
        #    tmpi = i * np.ones((len(np.arange(b, e)), ncomp))

        #    write_data[j:j+len(tmp), 0] = tmp
        #    write_data[j:j+len(tmp), 1] = tmpl
        #    write_data[j:j+len(tmp), 2:] = tmpi
        #    j += len(tmp)

#    def create_traj(self, top_n=None):
#        r"""
#        Create the customized trajectories for the individual mixture components
#        of the model. If `top_n` is None, a single trajectory is created
#        with all of `sel1` and a single `sel2` residue, with every contact
#        accounted for (ie. a single frame in the original trajectory may be
#        used multiple times due to multiple contacts formed with the `sel1`
#        residue of interest at that frame).
#
#        :param top_n: Number of frames desired for the individual trajectories
#                      (sorted in order of decreasing classification 
#                      probability).
#        :type top_n: int
#        """
#
#        if os.path.exists(self.fulltraj) and top_n is None:
#            raise FileExistsError(f'{self.fulltraj} exists, remove then rerun')
#
#        write_ag = self.ag1.atoms + self.ag2.residues[0].atoms
#        write_ag.atoms.write(self.topname)
#        if not os.path.exists(self.dataname):
#            self._create_data()
#
#        tmp = np.load(self.dataname, mmap_mode='r')
#        u = mda.Universe(f'{self.utop}', f'{self.utraj}')
#        ag1 = u.atoms[self.ag1.indices]
#        ag2 = u.atoms[self.ag2.indices]
#        if top_n is not None:
#            sortinds = [tmp[:, i+2].argsort()[::-1][:top_n] for i in
#                        range(self.gibbs.processed_results.ncomp)]
#            for k in range(self.gibbs.processed_results.ncomp):
#                swf = tmp[sortinds[k], 0].astype(int)
#                swl = tmp[sortinds[k], 1].astype(int)
#                with mda.Writer(f'basicrta-{self.cutoff}/{self.gibbs.residue}/'
#                                f'chol_traj_comp{k}_top{top_n}.xtc',
#                                len(write_ag.atoms)) as W:
#                    for i, ts in tqdm(enumerate(u.trajectory[swf]),
#                                      total=len(swf),
#                                      desc=f'writing component {k}'):
#                        W.write(ag1 + ag2.select_atoms(f'resid {swl[i]}'))
#
#        else:
#            with mda.Writer(self.fulltraj, len(write_ag.atoms)) as W:
#                for i, ts in tqdm(enumerate(u.trajectory[tmp[:, 0].
#                                  astype(int)]), total=len(tmp),
#                                  desc='writing trajectory'):
#                    W.write(ag1 + ag2.select_atoms(f'resid {int(tmp[i, 1])}'))

def get_parser():
    import argparse
    parser = argparse.ArgumentParser(description="""map kinetics from clustered
                                     results onto trajectory, create weighted
                                     densities if flag is used""")
    required = parser.add_argument_group('required arguments')
    required.add_argument("--gibbs", type=str, required=True, help="""gibbs pickle
                        file to use for creating kinetic trajectories and
                        densities""")
    required.add_argument("--contacts", type=str, required=True, help="""contacts
                        file used in creation of the gibbs sampler data""")
    parser.add_argument("--top_n", type=int, nargs='?', help="""use the `top_n`
                        most likely frames to create trajectory or densities""")
    parser.add_argument("--step", type=int, nargs='?', default=1, help="""write
                        out frame if frame%%step=0""")
    parser.add_argument("--wdensity", action='store_true', help="""create
                        weighted densities""")
    # this is to make the cli work, should be just a temporary solution
    parser.add_argument('kinetics', nargs='?', help=argparse.SUPPRESS)
    return parser

def main():
    from basicrta.gibbs import Gibbs
    parser = get_parser()
    args = parser.parse_args()

    g = Gibbs().load(args.gibbs)
    mk = MapKinetics(g, args.contacts)
    mk.create_traj(top_n=args.top_n)
    if args.wdensity:
        mk.weighted_densities(step=args.step, top_n=args.top_n)

if __name__ == "__main__":
    exit(main())
