from tqdm import tqdm
from MDAnalysis.lib import distances
from multiprocessing import Pool, Lock
from basicrta import istarmap
import numpy as np
import multiprocessing
import collections
import MDAnalysis as mda
import pickle
import glob
import os
os.environ['MKL_NUM_THREADS'] = '1'


class MapContacts(object):
    """This class is used to create the map of contacts between two groups of
    atoms. A single cutoff is used to define a contact between the two groups,
    where if any atomic distance between the two groups is less than the cutoff,
    a contact is considered formed.
    
    :param u: Universe containing the topology and trajectory for which the
              contacts will be computed.
    :type u: `MDAnalysis Universe`
    :param ag1: Primary AtomGroup for which contacts will be computed, typically a
                protein.
    :type ag1: MDAnalysis AtomGroup
    :param ag2: Secondary AtomGroup which forms contacts with `ag1`, typically
                lipids, ions, or other small molecules. Each residue of `ag2` 
                must have the same number of atoms.
    :type ag2: MDAnalysis AtomGroup
    :param nproc: Number of processes to use in computing contacts (default is
                  1).
    :type nproc: int, optional
    :param frames: List of frames to use in computing contacts (default is
                   None, meaning all frames are used).
    :type frames: list or np.array, optional
    :param max_cutoff: Maximum cutoff to use in computing contacts. A primary 
                   contact map is created upon which multiple cutoffs can be
                   imposed, i.e. in the case where a proper cutoff is being
                   determined. This can typically be left at its default value,
                   unless a greater value is needed (default is 10.0).
    :type max_cutoff: float, optional
    :param nslices: Number of slices to break the trajectory into for
                    processing. If device memory is limited, try increasing
                    `nslices` (default is 100).
    :type nslices: int, optional
    """

    def __init__(self, u, ag1, ag2, nproc=1, frames=None, 
                 max_cutoff=10.0, nslices=100):
        self.u, self.nproc = u, nproc
        self.ag1, self.ag2 = ag1, ag2
        self.max_cutoff = max_cutoff
        self.frames, self.nslices = frames, nslices
        self.contacts_filename = f"contacts_max{self.max_cutoff}.pkl"

    def run(self):
        """Run contact analysis and save to `contacts_max{max_cutoff}.pkl`
        """
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
                                   'max_cutoff': self.max_cutoff})

        contact_map = np.memmap('.tmpmap', mode='w+',
                                shape=(mapsize, 5), dtype=dtype)
        for i in range(self.nslices):
            contact_map[bounds[i]:bounds[i+1]] = np.genfromtxt(f'.contacts_'
                                                               f'{i:04}',
                                                               delimiter=',')
            contact_map.flush()

        contact_map.dump(self.contacts_filename, protocol=5)
        os.remove('.tmpmap')
        cfiles = glob.glob('.contacts*')
        [os.remove(f) for f in cfiles]
        print(f'\nSaved contacts as "{self.contacts_filename}"')

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
                                              max_cutoff=self.max_cutoff)
                pairlist = [(self.ag1.resids[b[0][i, 0]],
                             self.ag2.resids[b[0][i, 1]]) for i in
                            range(len(b[0]))]
                pairdir = collections.Counter(a for a in pairlist)
                lsum = 0
                for j in pairdir:
                    temp = pairdir[j]
                    dset.append([int(ts.frame), int(j[0]), int(j[1]),
                                 float(min(b[1][lsum:lsum+temp])),
                                 float(np.round(ts.time, dec)/1000)])  # convert to ns
                    lsum += temp
                [f.write(f"{line}".strip('[]') + "\n") for line in dset]
                data_len += len(dset)
            f.flush()
        return data_len


class ProcessContacts(object):
    """The :class:`ProcessProtein` class takes the primary contact map
    (i.e. `contacts_max10.0.pkl`) and collects contacts based on a prescribed 
    cutoff. 

    :param cutoff: Collect all contacts between `ag1` and `ag2` within this
                   value.
    :type cutoff: float
    :param map_name: Name of primary contact map. The default produced by
                     MapContacts is `contacts_max10.0.pkl`.
    :type map_name: str
    :param nproc: Number of processes to use in collecting contacts (default is
                  1). 
    :type nproc: int, optional
    """

    def __init__(self, cutoff, map_name, nproc=1):
        self.nproc = nproc
        self.map_name = map_name
        self.cutoff = cutoff

    def run(self):
        """Process contacts using the prescribed cutoff and write to
           contacts-{cutoff}.pkl
        """
        if os.path.exists(self.map_name):
            with open(self.map_name, 'r+b') as f:
                memmap = pickle.load(f)
            # memmap = np.load(self.map_name, mmap_mode='r')
            dtype = memmap.dtype
            new_metadata = dtype.metadata.copy()
            new_metadata['cutoff'] = self.cutoff
            new_dtype = np.dtype(np.float64, metadata=new_metadata)

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
        contact_map = np.memmap('.tmpmap', mode='w+',
                                shape=(mapsize, 4), dtype=new_dtype)

        for i in range(len(lresids)):
            contact_map[bounds[i]:bounds[i+1]] = np.load(f'.contacts_{i:04}.'
                                                         f'npy')
            contact_map.flush()

        contact_map.dump(f'contacts_{self.cutoff}.pkl', protocol=5)
        # os.remove('.tmpmap')
        # cfiles = glob.glob('.contacts*')
        # [os.remove(f) for f in cfiles]
        print(f'\nSaved contacts to "contacts_{self.cutoff}.pkl"')

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


class CombineContacts(object):
    """Class to combine contact timeseries from multiple repeat runs.
    
    This class enables pooling data from multiple trajectory repeats and
    calculating posteriors from all data together, rather than analyzing
    each run separately.
    
    :param contact_files: List of contact pickle files to combine
    :type contact_files: list of str
    :param output_name: Name for the combined output file (default: 'combined_contacts.pkl')
    :type output_name: str, optional
    :param validate_compatibility: Whether to validate that files are compatible (default: True)
    :type validate_compatibility: bool, optional
    """
    
    def __init__(self, contact_files, output_name='combined_contacts.pkl', 
                 validate_compatibility=True):
        self.contact_files = contact_files
        self.output_name = output_name
        self.validate_compatibility = validate_compatibility
        
        if len(contact_files) < 2:
            raise ValueError("At least 2 contact files are required for combining")
            
    def _load_contact_file(self, filename):
        """Load a contact pickle file and return data and metadata."""
        if not os.path.exists(filename):
            raise FileNotFoundError(f"Contact file not found: {filename}")
            
        with open(filename, 'rb') as f:
            contacts = pickle.load(f)
            
        metadata = contacts.dtype.metadata
        return contacts, metadata
        
    def _validate_compatibility(self, metadatas):
        """Validate that contact files are compatible for combining."""
        reference = metadatas[0]
        
        # Check that all files have the same atom groups
        for i, meta in enumerate(metadatas[1:], 1):
            # Compare cutoff
            if meta['cutoff'] != reference['cutoff']:
                raise ValueError(f"Incompatible cutoffs: file 0 has {reference['cutoff']}, "
                               f"file {i} has {meta['cutoff']}")
                               
            # Compare atom group selections by checking if resids match
            ref_ag1_resids = set(reference['ag1'].residues.resids)
            ref_ag2_resids = set(reference['ag2'].residues.resids)
            meta_ag1_resids = set(meta['ag1'].residues.resids)
            meta_ag2_resids = set(meta['ag2'].residues.resids)
            
            if ref_ag1_resids != meta_ag1_resids:
                raise ValueError(f"Incompatible ag1 residues between file 0 and file {i}")
            if ref_ag2_resids != meta_ag2_resids:
                raise ValueError(f"Incompatible ag2 residues between file 0 and file {i}")
                
        # Check timesteps and warn if different
        timesteps = [meta['ts'] for meta in metadatas]
        if not all(abs(ts - timesteps[0]) < 1e-6 for ts in timesteps):
            print("WARNING: Different timesteps detected across runs:")
            for i, (filename, ts) in enumerate(zip(self.contact_files, timesteps)):
                print(f"  File {i} ({filename}): dt = {ts} ns")
            print("This may affect residence time estimates, especially for fast events.")
            
    def run(self):
        """Combine contact files and save the result."""
        print(f"Combining {len(self.contact_files)} contact files...")
        
        all_contacts = []
        all_metadatas = []
        
        # Load all contact files
        for i, filename in enumerate(self.contact_files):
            print(f"Loading file {i+1}/{len(self.contact_files)}: {filename}")
            contacts, metadata = self._load_contact_file(filename)
            all_contacts.append(contacts)
            all_metadatas.append(metadata)
            
        # Validate compatibility if requested
        if self.validate_compatibility:
            print("Validating file compatibility...")
            self._validate_compatibility(all_metadatas)
            
        # Combine contact data
        print("Combining contact data...")
        
        # Calculate total size and create combined array
        total_size = sum(len(contacts) for contacts in all_contacts)
        reference_metadata = all_metadatas[0].copy()
        
        # Extend metadata to include trajectory source information
        reference_metadata['source_files'] = self.contact_files
        reference_metadata['n_trajectories'] = len(self.contact_files)
        
        # Determine number of columns (5 for raw contacts, 4 for processed)
        n_cols = all_contacts[0].shape[1]
        
        # Create dtype with extended metadata
        combined_dtype = np.dtype(np.float64, metadata=reference_metadata)
        
        # Add trajectory source column (will be last column)
        combined_contacts = np.zeros((total_size, n_cols + 1), dtype=np.float64)
        
        # Combine data and add trajectory source information
        offset = 0
        for traj_idx, contacts in enumerate(all_contacts):
            n_contacts = len(contacts)
            # Copy original contact data
            combined_contacts[offset:offset+n_contacts, :n_cols] = contacts[:]
            # Add trajectory source index
            combined_contacts[offset:offset+n_contacts, n_cols] = traj_idx
            offset += n_contacts
            
        # Create final memmap with proper dtype
        final_contacts = combined_contacts.view(combined_dtype)
        
        # Save combined contacts
        print(f"Saving combined contacts to {self.output_name}...")
        final_contacts.dump(self.output_name, protocol=5)
        
        print(f"Successfully combined {len(self.contact_files)} files into {self.output_name}")
        print(f"Total contacts: {total_size}")
        print(f"Added trajectory source column (index {n_cols}) for kinetic clustering support")
        
        return self.output_name


if __name__ == '__main__':
    """DOCSSS
    """
    import argparse
    parser = argparse.ArgumentParser(description="Create the primary contact \
                                     map and collect contacts based on the \
                                     desired cutoff distance")
    parser.add_argument('--top', type=str, help="Topology")
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
    
    mc = MapContacts(u, ag1, ag2, nproc=nproc, nslices=nslices)
    mapname = mc.contacts_filename
    if not os.path.exists(mapname):
        print(f"running MapContacts to generate {mapname}")
        mc.run()
    else:
        print(f"using existing {mapname}")

    ProcessContacts(cutoff, mapname, nproc=nproc).run()

