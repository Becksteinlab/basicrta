import warnings
import pytest 
import pickle
import MDAnalysis as mda
import numpy as np
import os
import sys

@pytest.mark.filterwarnings("ignore::UserWarning")
def test_mapcontacts():
    from basicrta.tests.datafiles import PDB, XTC
    from basicrta.contacts import MapContacts

    u = mda.Universe(PDB, XTC)
    P88 = u.select_atoms('resname PRO and resid 88')
    chol = u.select_atoms('resname CHOL and resid 309')

    MapContacts(u, P88, chol, nslices=1).run()

def test_contacts():
    with open('basicrta/tests/contacts.pkl', 'rb') as c:
        contacts = pickle.load(c)
    
    filtered_contacts = contacts[contacts[:,3] <= 7]
    assert len(filtered_contacts) == 5
    assert (filtered_contacts[:,0] == [96,97,98,99,100]).all() 

@pytest.mark.filterwarnings("ignore::UserWarning")
def test_max_cutoff():
    from basicrta.tests.datafiles import PDB, XTC
    from basicrta.contacts import MapContacts

    u = mda.Universe(PDB, XTC)
    P88 = u.select_atoms('resname PRO and resid 88')
    chol = u.select_atoms('resname CHOL and resid 309')

    MapContacts(u, P88, chol, nslices=1, max_cutoff=12.0).run()
    
    with open('basicrta/tests/contacts.pkl', 'rb') as c:
        contacts = pickle.load(c)
    
    assert len(contacts) == 30
    assert (contacts[:, 0] == np.delete(np.arange(69,101), [4,5])).all()

def test_contact_metadata():
    with open('basicrta/tests/contacts.pkl', 'rb') as c:
        contacts = pickle.load(c)

    assert list(contacts.dtype.metadata) == ['top', 'traj', 'ag1', 'ag2', 'ts',
                                             'max_cutoff']

def test_processcontacts():
    from basicrta.contacts import ProcessContacts
    ProcessContacts(7.0).run()

def test_processed_contacts():
    with open('basicrta/tests/contacts_7.0.pkl', 'rb') as f:
       contacts = pickle.load(f)

    assert (contacts == [88, 309, 9.6, 0.5]).all()

def test_processed_contact_metadata():
    with open('basicrta/tests/contacts_7.0.pkl', 'rb') as c:
        contacts = pickle.load(c)

    assert list(contacts.dtype.metadata) == ['top', 'traj', 'ag1', 'ag2', 'ts',
                                             'max_cutoff', 'cutoff']
def test_tearDown():
    import os
    os.remove('contacts.pkl')
    os.remove('contacts_7.0.pkl')
    os.remove('.tmpmap')

