"""
Tests for combining contact timeseries from multiple repeat runs.
"""

import os
import pytest
import numpy as np
import pickle

from basicrta.contacts import CombineContacts


class MockAtomGroup:
    """Mock atom group that can be pickled."""
    def __init__(self, resids):
        self.residues = MockResidues(resids)
        

class MockResidues:
    """Mock residues that can be pickled."""
    def __init__(self, resids):
        self.resids = np.array(resids)





@pytest.fixture
def create_mock_contacts():
    """Factory fixture for creating mock contact files."""
    def _create_mock_contacts(tmp_path, filename, n_contacts=100, cutoff=7.0, 
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
        rng = np.random.default_rng(seed=42)  # For reproducible tests
        contacts = np.zeros((n_contacts, 4), dtype=dtype)
        contacts[:, 0] = rng.choice([1, 3, 5], n_contacts)  # protein resids
        contacts[:, 1] = rng.choice([100, 101, 102, 201, 202], n_contacts)  # lipid resids  
        contacts[:, 2] = rng.uniform(0, 100, n_contacts)  # start times
        
        # Sample from a 2-term hyperexponential distribution
        # for the waiting times t:
        # t ~ 0.7 * exp(-10.0 * t) + 0.3 * exp(-0.5 * t)
        u = rng.uniform(0, 1, n_contacts)
        component_choice = rng.choice([0, 1], n_contacts, p=[0.7, 0.3])
        residence_times = np.where(
            component_choice == 0,
            rng.exponential(1/10.0, n_contacts),  # First component: rate 10.0
            rng.exponential(1/0.5, n_contacts)    # Second component: rate 0.5
        )
        contacts[:, 3] = residence_times
        
        # Save to file
        filepath = tmp_path / filename
        with open(filepath, 'wb') as f:
            pickle.dump(contacts, f, protocol=5)
            
        return contacts, metadata
    
    return _create_mock_contacts


class TestCombineContacts:
    """Test class for CombineContacts functionality."""
        
    def test_combine_contacts_basic(self, tmp_path, create_mock_contacts):
        """Test basic contact combination functionality."""
        # Create two mock contact files
        contacts1, meta1 = create_mock_contacts(tmp_path, "contacts1.pkl", n_contacts=50)
        contacts2, meta2 = create_mock_contacts(tmp_path, "contacts2.pkl", n_contacts=75, 
                                                    traj_name="test2.xtc")
        
        # Combine them
        combiner = CombineContacts(
            contact_files=[str(tmp_path / "contacts1.pkl"), str(tmp_path / "contacts2.pkl")],
            output_name=str(tmp_path / "combined.pkl")
        )
        
        output_file = combiner.run()
        
        # Verify output
        assert output_file == str(tmp_path / "combined.pkl")
        assert os.path.exists(tmp_path / "combined.pkl")
        
        # Load and verify combined data
        with open(tmp_path / "combined.pkl", 'rb') as f:
            combined = pickle.load(f)
            
        # Check shape - should have original 4 cols + 1 trajectory source col  
        assert combined.shape == (125, 5)
        
        # Check metadata
        metadata = combined.dtype.metadata
        assert metadata['source_files'] == [str(tmp_path / "contacts1.pkl"), str(tmp_path / "contacts2.pkl")]
        assert metadata['n_trajectories'] == 2
        assert metadata['cutoff'] == 7.0
        
        # Check trajectory source column (last column)
        traj_sources = combined[:, 4]
        assert np.all(traj_sources[:50] == 0)  # First 50 from file 0
        assert np.all(traj_sources[50:] == 1)  # Next 75 from file 1
        
    def test_incompatible_cutoffs(self, tmp_path, create_mock_contacts):
        """Test that incompatible cutoffs raise an error."""
        create_mock_contacts(tmp_path, "contacts1.pkl", cutoff=7.0)
        create_mock_contacts(tmp_path, "contacts2.pkl", cutoff=8.0)  # Different cutoff
        
        combiner = CombineContacts(
            contact_files=[str(tmp_path / "contacts1.pkl"), str(tmp_path / "contacts2.pkl")]
        )
        
        with pytest.raises(ValueError, match="Incompatible cutoffs"):
            combiner.run()
            
    def test_incompatible_atom_groups(self, tmp_path, create_mock_contacts):
        """Test that incompatible atom groups raise an error."""
        # Create first file with standard residues
        contacts1, _ = create_mock_contacts(tmp_path, "contacts1.pkl")
        
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
        
        with open(tmp_path / "contacts2.pkl", 'wb') as f:
            pickle.dump(contacts, f, protocol=5)
        
        combiner = CombineContacts(
            contact_files=[str(tmp_path / "contacts1.pkl"), str(tmp_path / "contacts2.pkl")]
        )
        
        with pytest.raises(ValueError, match="Incompatible ag1 residues"):
            combiner.run()
            
    def test_different_timesteps_warning(self, tmp_path, create_mock_contacts, capsys):
        """Test that different timesteps produce a warning."""
        create_mock_contacts(tmp_path, "contacts1.pkl", ts=0.1)
        create_mock_contacts(tmp_path, "contacts2.pkl", ts=0.2)  # Different timestep
        
        combiner = CombineContacts(
            contact_files=[str(tmp_path / "contacts1.pkl"), str(tmp_path / "contacts2.pkl")]
        )
        
        combiner.run()
        
        # Check that warning was printed
        captured = capsys.readouterr()
        assert "WARNING: Different timesteps detected" in captured.out
        
    def test_minimum_files_required(self, tmp_path, create_mock_contacts):
        """Test that at least 2 files are required."""
        create_mock_contacts(tmp_path, "contacts1.pkl")
        
        with pytest.raises(ValueError, match="At least 2 contact files are required"):
            CombineContacts(contact_files=[str(tmp_path / "contacts1.pkl")])
            
    def test_missing_file(self, tmp_path, create_mock_contacts):
        """Test handling of missing contact files."""
        create_mock_contacts(tmp_path, "contacts1.pkl")
        
        combiner = CombineContacts(
            contact_files=[str(tmp_path / "contacts1.pkl"), str(tmp_path / "nonexistent.pkl")]
        )
        
        with pytest.raises(FileNotFoundError, match="Contact file not found"):
            combiner.run()
            
    def test_skip_validation(self, tmp_path, create_mock_contacts):
        """Test skipping compatibility validation."""
        create_mock_contacts(tmp_path, "contacts1.pkl", cutoff=7.0)
        create_mock_contacts(tmp_path, "contacts2.pkl", cutoff=8.0)  # Different cutoff
        
        combiner = CombineContacts(
            contact_files=[str(tmp_path / "contacts1.pkl"), str(tmp_path / "contacts2.pkl")],
            validate_compatibility=False
        )
        
        # Should not raise error when validation is skipped
        output_file = combiner.run()
        assert os.path.exists(output_file)
        
    def test_combined_contacts_detection(self, tmp_path, create_mock_contacts):
        """Test that combined contact files are properly detected."""
        # Create and combine contacts
        create_mock_contacts(tmp_path, "contacts1.pkl", n_contacts=30)
        create_mock_contacts(tmp_path, "contacts2.pkl", n_contacts=40, traj_name="test2.xtc")
        
        combiner = CombineContacts(
            contact_files=[str(tmp_path / "contacts1.pkl"), str(tmp_path / "contacts2.pkl")],
            output_name=str(tmp_path / "combined.pkl")
        )
        
        combiner.run()
        
        # Load combined file and check metadata
        with open(tmp_path / "combined.pkl", 'rb') as f:
            combined = pickle.load(f)
            
        metadata = combined.dtype.metadata
        assert 'n_trajectories' in metadata
        assert metadata['n_trajectories'] == 2
        assert 'source_files' in metadata
        assert len(metadata['source_files']) == 2

    def test_gibbs_sampler_integration(self, tmp_path, create_mock_contacts):
        """Test that Gibbs sampler works with combined timeseries."""
        from basicrta.gibbs import Gibbs
        
        # Create mock contact files with more realistic data for Gibbs sampling
        create_mock_contacts(tmp_path, "contacts1.pkl", n_contacts=100, cutoff=7.0)
        create_mock_contacts(tmp_path, "contacts2.pkl", n_contacts=100, cutoff=7.0, traj_name="test2.xtc")
        
        # Combine the contacts
        combiner = CombineContacts(
            contact_files=[str(tmp_path / "contacts1.pkl"), str(tmp_path / "contacts2.pkl")],
            output_name=str(tmp_path / "combined.pkl")
        )
        
        combined_file = combiner.run()
        
        # Load combined contacts
        with open(combined_file, 'rb') as f:
            combined_contacts = pickle.load(f)
        
        # Extract residence times for Gibbs sampling
        residence_times = combined_contacts[:, 3]  # 4th column is residence time
        
        # Create a temporary directory for Gibbs output
        output_dir = tmp_path / "basicrta-7.0" / "test_residue"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Run Gibbs sampler with minimal settings
        gibbs = Gibbs(
            times=residence_times,
            residue='test_residue',
            ncomp=2,  # Use 2 components for stability
            niter=1000,  # 1000 steps as requested
            burnin=5,   # 5 burnin steps as requested
            cutoff=7.0  # Set cutoff for directory creation
        )
        
        # Run the sampler (using Python API as requested)
        gibbs.run()
        
        # Verify that results were generated
        assert hasattr(gibbs, 'mcrates')
        assert hasattr(gibbs, 'mcweights')
        assert gibbs.mcrates is not None
        assert gibbs.mcweights is not None
        
        # Check that we have the expected number of iterations stored
        # (samples are saved every g steps, so 1000/50 = 20 samples)
        assert len(gibbs.mcrates) == 20  # 1000 iterations / 50 g steps = 20 samples
        assert len(gibbs.mcweights) == 20
        
        # Verify shapes are correct
        assert gibbs.mcrates.shape == (20, 2)  # 20 samples, 2 components
        assert gibbs.mcweights.shape == (20, 2)  # 20 samples, 2 components