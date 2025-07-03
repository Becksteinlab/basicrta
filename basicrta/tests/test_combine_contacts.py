"""
Tests for combining contact timeseries from multiple repeat runs.
"""

import os
import pytest
import numpy as np
import pickle
import tempfile
import shutil

from basicrta.contacts import CombineContacts


class MockAtomGroup:
    """Mock atom group that can be pickled."""
    def __init__(self, resids):
        self.residues = MockResidues(resids)
        

class MockResidues:
    """Mock residues that can be pickled."""
    def __init__(self, resids):
        self.resids = np.array(resids)


class TestCombineContacts:
    """Test class for CombineContacts functionality."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.original_dir = os.getcwd()
        os.chdir(self.temp_dir)
        
    def teardown_method(self):
        """Clean up after tests."""
        os.chdir(self.original_dir)
        shutil.rmtree(self.temp_dir)
        
    def create_mock_contacts(self, filename, n_contacts=100, cutoff=7.0, 
                           ts=0.1, traj_name="test.xtc", top_name="test.pdb"):
        """Create a mock contact file for testing."""
        # Create mock atom groups that can be pickled
        mock_ag1 = MockAtomGroup([1, 2, 3, 4, 5])
        mock_ag2 = MockAtomGroup([100, 101, 102])
        
        # Create metadata
        metadata = {
            'top': top_name,
            'traj': traj_name,
            'ag1': mock_ag1,
            'ag2': mock_ag2, 
            'ts': ts,
            'cutoff': cutoff
        }
        
        # Create dtype with metadata
        dtype = np.dtype(np.float64, metadata=metadata)
        
        # Create contact data (4 columns for processed contacts)
        # [protein_resid, lipid_resid, start_time, residence_time]
        np.random.seed(42)  # For reproducible tests
        contacts = np.zeros((n_contacts, 4), dtype=dtype)
        contacts[:, 0] = np.random.choice([1, 2, 3, 4, 5], n_contacts)  # protein resids
        contacts[:, 1] = np.random.choice([100, 101, 102], n_contacts)  # lipid resids  
        contacts[:, 2] = np.random.uniform(0, 100, n_contacts)  # start times
        contacts[:, 3] = np.random.exponential(1.0, n_contacts)  # residence times
        
        # Save to file
        with open(filename, 'wb') as f:
            pickle.dump(contacts, f, protocol=5)
            
        return contacts, metadata
        
    def test_combine_contacts_basic(self):
        """Test basic contact combination functionality."""
        # Create two mock contact files
        contacts1, meta1 = self.create_mock_contacts("contacts1.pkl", n_contacts=50)
        contacts2, meta2 = self.create_mock_contacts("contacts2.pkl", n_contacts=75, 
                                                    traj_name="test2.xtc")
        
        # Combine them
        combiner = CombineContacts(
            contact_files=["contacts1.pkl", "contacts2.pkl"],
            output_name="combined.pkl"
        )
        
        output_file = combiner.run()
        
        # Verify output
        assert output_file == "combined.pkl"
        assert os.path.exists("combined.pkl")
        
        # Load and verify combined data
        with open("combined.pkl", 'rb') as f:
            combined = pickle.load(f)
            
        # Check shape - should have original 4 cols + 1 trajectory source col  
        assert combined.shape == (125, 5)
        
        # Check metadata
        metadata = combined.dtype.metadata
        assert metadata['source_files'] == ["contacts1.pkl", "contacts2.pkl"]
        assert metadata['n_trajectories'] == 2
        assert metadata['cutoff'] == 7.0
        
        # Check trajectory source column (last column)
        traj_sources = combined[:, 4]
        assert np.all(traj_sources[:50] == 0)  # First 50 from file 0
        assert np.all(traj_sources[50:] == 1)  # Next 75 from file 1
        
    def test_incompatible_cutoffs(self):
        """Test that incompatible cutoffs raise an error."""
        self.create_mock_contacts("contacts1.pkl", cutoff=7.0)
        self.create_mock_contacts("contacts2.pkl", cutoff=8.0)  # Different cutoff
        
        combiner = CombineContacts(
            contact_files=["contacts1.pkl", "contacts2.pkl"]
        )
        
        with pytest.raises(ValueError, match="Incompatible cutoffs"):
            combiner.run()
            
    def test_incompatible_atom_groups(self):
        """Test that incompatible atom groups raise an error."""
        # Create first file with standard residues
        contacts1, _ = self.create_mock_contacts("contacts1.pkl")
        
        # Create second file with different protein residues
        mock_ag1 = MockAtomGroup([10, 20, 30])  # Different resids
        mock_ag2 = MockAtomGroup([100, 101, 102])
        
        metadata = {
            'top': "test2.pdb",
            'traj': "test2.xtc", 
            'ag1': mock_ag1,
            'ag2': mock_ag2,
            'ts': 0.1,
            'cutoff': 7.0
        }
        
        dtype = np.dtype(np.float64, metadata=metadata)
        contacts = np.zeros((50, 4), dtype=dtype)
        
        with open("contacts2.pkl", 'wb') as f:
            pickle.dump(contacts, f, protocol=5)
        
        combiner = CombineContacts(
            contact_files=["contacts1.pkl", "contacts2.pkl"]
        )
        
        with pytest.raises(ValueError, match="Incompatible ag1 residues"):
            combiner.run()
            
    def test_different_timesteps_warning(self, capsys):
        """Test that different timesteps produce a warning."""
        self.create_mock_contacts("contacts1.pkl", ts=0.1)
        self.create_mock_contacts("contacts2.pkl", ts=0.2)  # Different timestep
        
        combiner = CombineContacts(
            contact_files=["contacts1.pkl", "contacts2.pkl"]
        )
        
        combiner.run()
        
        # Check that warning was printed
        captured = capsys.readouterr()
        assert "WARNING: Different timesteps detected" in captured.out
        
    def test_minimum_files_required(self):
        """Test that at least 2 files are required."""
        self.create_mock_contacts("contacts1.pkl")
        
        with pytest.raises(ValueError, match="At least 2 contact files are required"):
            CombineContacts(contact_files=["contacts1.pkl"])
            
    def test_missing_file(self):
        """Test handling of missing contact files."""
        self.create_mock_contacts("contacts1.pkl")
        
        combiner = CombineContacts(
            contact_files=["contacts1.pkl", "nonexistent.pkl"]
        )
        
        with pytest.raises(FileNotFoundError, match="Contact file not found"):
            combiner.run()
            
    def test_skip_validation(self):
        """Test skipping compatibility validation."""
        self.create_mock_contacts("contacts1.pkl", cutoff=7.0)
        self.create_mock_contacts("contacts2.pkl", cutoff=8.0)  # Different cutoff
        
        combiner = CombineContacts(
            contact_files=["contacts1.pkl", "contacts2.pkl"],
            validate_compatibility=False
        )
        
        # Should not raise error when validation is skipped
        output_file = combiner.run()
        assert os.path.exists(output_file)