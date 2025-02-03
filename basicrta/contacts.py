#!/usr/bin/env python

import os
os.environ['MKL_NUM_THREADS'] = '1'
from tqdm import tqdm
import numpy as np
import multiprocessing
from MDAnalysis.lib import distances
import collections
from multiprocessing import Pool, Lock
import MDAnalysis as mda
import pickle
import glob
from basicrta import istarmap


class MapContacts(object):
    """
    This class is used to create the map of contacts between two groups of
    atoms. A single cutoff is used to define a contact between the two groups,
    where if any atomic distance between the two groups is less than the cutoff,
    a contact is considered formed.
    """

    def __init__(self, u, ag1, ag2, nproc=1, frames=None, cutoff=10.0,
                 nslices=100):
        self.u, self.nproc = u, nproc
        self.ag1, self.ag2 = ag1, ag2
        self.cutoff, self.frames, self.nslices = cutoff, frames, nslices

    def run(self):
        if self.frames is not None:
            sliced_frames = np.array_split(self.frames, self.nslices)
        else:
            sliced_frames = np.array_split(np.arange(len(self.u.trajectory)),
                                           self.nslices)

        input_list = [[i, self.u.trajectory[aslice]] for
                      i, aslice in enumerate(sliced_frames)]

        lens = []
        with (Pool(self.nproc, initializer=tqdm.set_lock, initargs=(Lock(),))
              as p):
            for alen in tqdm(p.istarmap(self._run_contacts, input_list),
                          total=self.nslices, position=0,
                          desc='overall progress'):
                lens.append(alen)
        lens = np.array(lens)
        mapsize = sum(lens)
        bounds = np.concatenate([[0], np.cumsum(lens)])
        dtype = np.dtype(np.float64,
                         metadata={'top': self.u.filename,
                                   'traj': self.u.trajectory.filename,
                                   'ag1': self.ag1, 'ag2': self.ag2,
                                   'ts': self.u.trajectory.dt/1000,
                                   'cutoff': self.cutoff})

        contact_map = np.memmap('.tmpmap', mode='w+',
                                shape=(mapsize, 5), dtype=dtype)
        for i in range(self.nslices):
            contact_map[bounds[i]:bounds[i+1]] = np.genfromtxt(f'.contacts_'
                                                               f'{i:04}',
                                                               delimiter=',')
            contact_map.flush()

        contact_map.dump('contacts.pkl', protocol=5)
        os.remove('.tmpmap')
        cfiles = glob.glob('.contacts*')
        [os.remove(f) for f in cfiles]
        print('\nSaved contacts as "contacts.pkl"')

    def _run_contacts(self, i, sliced_traj):
        from basicrta.util import get_dec

        try:
            proc = int(multiprocessing.current_process().name.split('-')[-1])
        except ValueError:
            proc = 1

        with open(f'.contacts_{i:04}', 'w+') as f:
            dec = get_dec(self.u.trajectory.ts.dt/1000)  # convert to ns
            text = f'slice {i+1} of {self.nslices}'
            data_len = 0
            for ts in tqdm(sliced_traj, desc=text, position=proc,
                           total=len(sliced_traj), leave=False):
                dset = []
                b = distances.capped_distance(self.ag1.positions,
                                              self.ag2.positions,
                                              max_cutoff=self.cutoff)
                pairlist = [(self.ag1.resids[b[0][i, 0]],
                             self.ag2.resids[b[0][i, 1]]) for i in
                            range(len(b[0]))]
                pairdir = collections.Counter(a for a in pairlist)
                lsum = 0
                for j in pairdir:
                    temp = pairdir[j]
                    dset.append([ts.frame, j[0], j[1],
                                 min(b[1][lsum:lsum+temp]),
                                 np.round(ts.time, dec)/1000])  # convert to ns
                    lsum += temp
                [f.write(f"{line}".strip('[]') + "\n") for line in dset]
                data_len += len(dset)
            f.flush()
        return data_len


class ProcessContacts(object):
    def __init__(self, cutoff, nproc, map_name='contacts.pkl'):
        self.nproc = nproc
        self.map_name = map_name
        self.cutoff = cutoff

    def run(self):
        from basicrta.util import siground

        if os.path.exists(self.map_name):
            with open(self.map_name, 'r+b') as f:
                memmap = pickle.load(f)
            # memmap = np.load(self.map_name, mmap_mode='r')
            dtype = memmap.dtype
            memmap = memmap[memmap[:, -2] <= self.cutoff]
        else:
            raise FileNotFoundError(f'{self.map_name} not found. Specify the '
                                    'contacts file using the "map_name" '
                                    'argument')

        self.ts = dtype.metadata['ts']
        lresids = np.unique(memmap[:, 2])
        params = [[res, memmap[memmap[:, 2] == res], i] for i, res in
                  enumerate(lresids)]
        pool = Pool(self.nproc, initializer=tqdm.set_lock, initargs=(Lock(),))

        try:
            lens = pool.starmap(self._lipswap, params)
        except KeyboardInterrupt:
            pool.terminate()
        pool.close()

        bounds = np.concatenate([[0], np.cumsum(lens)]).astype(int)
        mapsize = sum(lens)
        contact_map = np.memmap(f'.tmpmap', mode='w+',
                                shape=(mapsize, 4), dtype=dtype)

        for i in range(len(lresids)):
            contact_map[bounds[i]:bounds[i+1]] = np.load(f'.contacts_{i:04}.'
                                                         f'npy')
            contact_map.flush()

        contact_map.dump(f'contacts_{self.cutoff}.pkl', protocol=5)
        # os.remove('.tmpmap')
        # cfiles = glob.glob('.contacts*')
        # [os.remove(f) for f in cfiles]
        print(f'\nSaved contacts to "contacts_{self.cutoff}.npy"')


    def _lipswap(self, lip, memarr, i):
        from basicrta.util import get_dec
        try:
            # proc = int(multiprocessing.current_process().name[-1])
            proc = int(multiprocessing.current_process().name.split('-')[-1])
        except ValueError:
            proc = 1

        presids = np.unique(memarr[:, 1])
        dset = []
        dec, ts = get_dec(self.ts), self.ts
        for pres in tqdm(presids, desc=f'lipID {int(lip)}', position=proc,
                         leave=False):
            stimes = np.round(memarr[:, -1][memarr[:, 1] == pres], dec)
            if len(stimes) == 0:
                continue
            stimes = np.concatenate([np.array([-1]), stimes,
                                     np.array([stimes[-1] + 1])])
            diff = np.round(stimes[1:] - stimes[:-1], dec)
            singles = stimes[
                np.where((diff[1:] > ts) & (diff[:-1] > ts))[0] + 1]
            diff[diff > ts] = 0
            inds = np.where(diff == 0)[0]
            sums = [sum(diff[inds[i]:inds[i + 1]]) for i in
                    range(len(inds) - 1)]
            clens = np.round(np.array(sums), dec)
            minds = np.where(clens != 0)[0]
            clens = clens[minds] + ts
            strt_times = stimes[inds[minds] + 1]

            [dset.append([pres, lip, time, ts]) for time in singles]
            [dset.append([pres, lip, time, clen]) for time, clen in
             zip(strt_times, clens)]
        np.save(f'.contacts_{i:04}', np.array(dset))
        return len(dset)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--top', type=str)
    parser.add_argument('--traj', type=str)
    parser.add_argument('--sel1', type=str)
    parser.add_argument('--sel2', type=str)
    parser.add_argument('--cutoff', type=float)
    parser.add_argument('--nproc', type=int, default=1)
    parser.add_argument('--nslices', type=int, default=100)
    args = parser.parse_args()

    u = mda.Universe(args.top, args.traj)
    cutoff, nproc, nslices = args.cutoff, args.nproc, args.nslices
    ag1 = u.select_atoms(args.sel1)
    ag2 = u.select_atoms(args.sel2)

    if not os.path.exists('contacts.pkl'):
        MapContacts(u, ag1, ag2, nproc=nproc, nslices=nslices).run()

    ProcessContacts(cutoff, nproc).run()
