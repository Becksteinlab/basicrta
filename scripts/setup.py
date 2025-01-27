import os 
import numpy as np
import MDAnalysis as mda
import pickle

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--contacts', type=str)
    args = parser.parse_args()

    dirname = args.contacts.strip('.pkl')
    cutoff = dirname.split('_')[1]

    with open(f'{args.contacts}', 'rb') as f:
        a = pickle.load(f)

    rg = a.dtype.metadata['ag1'].residues 
    ids = rg.resids
    names = rg.resnames
    names = np.array([mda.lib.util.convert_aa_code(name) for name in names])
    uniqs = np.unique(a[:, 0]).astype(int)

    inds = np.array([np.where(ids==val)[0][0] for val in uniqs])
    resids, resnames = ids[inds], names[inds]
    residues = np.array([f'{name}{resid}' for name, resid in zip(resnames, resids)])

    with open('residue_list.csv', 'w+') as r:
        for residue in residues:
            r.write(f'{residue},')

    if not os.path.exists(f'basicrta-{cutoff}'):
        os.mkdir(f'basicrta-{cutoff}')

    os.chdir(f'basicrta-{cutoff}')

    for residue in residues:
        if not os.path.exists(residue):
            os.mkdir(residue)

    
