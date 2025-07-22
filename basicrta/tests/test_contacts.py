import warnings
import pytest 
import pickle
import os
import sys
import MDAnalysis as mda
import numpy as np
from pathlib import Path
from basicrta.tests.datafiles import PDB, XTC
from basicrta.tests.utils import work_in
from basicrta.contacts import MapContacts

@pytest.fixture(scope="module")
def setup_mapcontacts(tmp_path_factory):
    u = mda.Universe(PDB, XTC)
    P88 = u.select_atoms('resname PRO and resid 88')
    chol = u.select_atoms('resname CHOL and resid 309')
    
    tmp_path_factory.mktemp("data")
    datadir = tmp_path_factory.mktemp("data")

    with work_in(datadir):
        MapContacts(u, P88, chol, nslices=1).run()
                                                                                             
    fn = datadir / "contacts_max10.0.pkl"
    assert fn.exists(), "Failed to locate {str(fn)}"
    return fn

@pytest.fixture
def setup_processcontacts(setup_mapcontacts, tmp_path_factory):
    from basicrta.contacts import ProcessContacts
    datadir = setup_mapcontacts.parents[0]

    with work_in(datadir):
        ProcessContacts(7.0, map_name='contacts_max10.0.pkl').run()
                                                                                             
    fn = datadir / "contacts_7.0.pkl"
    assert fn.exists(), "Failed to locate {str(fn)}"
    return fn

@pytest.mark.filterwarnings("ignore::UserWarning")
def test_mapcontacts(setup_mapcontacts):
    with open(setup_mapcontacts, 'rb') as c:
        contacts = pickle.load(c)

    print(contacts)
    filtered_contacts = contacts[contacts[:,3] <= 7]
    assert len(filtered_contacts) == 5
    assert (filtered_contacts[:,0] == [96,97,98,99,100]).all() 

@pytest.mark.filterwarnings("ignore::UserWarning")
def test_max_cutoff(tmp_path):
    with work_in(tmp_path):
        u = mda.Universe(PDB, XTC)
        P88 = u.select_atoms('resname PRO and resid 88')
        chol = u.select_atoms('resname CHOL and resid 309')
        MapContacts(u, P88, chol, nslices=1, max_cutoff=12.0).run()

        with open('contacts_max12.0.pkl', 'rb') as p:
            contacts = pickle.load(p)

        assert os.path.exists("contacts_max12.0.pkl"), "contacts_max12.0.pkl was not generated"
    assert len(contacts) == 30
    assert (contacts[:, 0] == np.delete(np.arange(69,101), [4,5])).all()

def test_contact_metadata(setup_mapcontacts):
    with open(setup_mapcontacts, 'rb') as c:
        contacts = pickle.load(c)

    assert list(contacts.dtype.metadata) == ['top', 'traj', 'ag1', 'ag2', 'ts',
                                             'max_cutoff']

def test_processcontacts(setup_processcontacts):
    with open(setup_processcontacts, 'rb') as f:
       contacts = pickle.load(f)

def test_processed_contacts(setup_processcontacts):
    with open(setup_processcontacts, 'rb') as f:
       contacts = pickle.load(f)

    assert (contacts == [88, 309, 9.6, 0.5]).all()

def test_processed_contact_metadata(setup_processcontacts):
    with open(setup_processcontacts, 'rb') as c:
        contacts = pickle.load(c)

    assert list(contacts.dtype.metadata) == ['top', 'traj', 'ag1', 'ag2', 'ts',
                                             'max_cutoff', 'cutoff']


