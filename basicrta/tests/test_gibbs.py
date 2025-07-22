"""
Tests for the Gibbs sampler with comprehensive coverage.
"""

import pytest
import numpy as np
import tempfile
import os
import pickle
from unittest.mock import Mock, patch
from basicrta.gibbs import Gibbs, ParallelGibbs
from basicrta.tests.utils import work_in
import MDAnalysis as mda


@pytest.fixture
def synthetic_timeseries():
    """
    Generate synthetic residence times from a bi-exponential distribution.
    
    Returns
    -------
    dict
        Dictionary containing test timeseries and expected parameters
    """
    rng = np.random.default_rng(seed=42)
    n_samples = 200
    
    # Create a bimodal distribution with known parameters
    # Fast component: rate ~2.0 (1/0.5), weight ~0.5
    times_short = rng.exponential(0.5, n_samples // 2)  
    # Slow component: rate ~0.2 (1/5.0), weight ~0.5  
    times_long = rng.exponential(5.0, n_samples // 2)
    
    test_times = np.concatenate([times_short, times_long])
    rng.shuffle(test_times)  # Mix them up
    
    return {
        'times': test_times,
        'expected_components': 2,
        'expected_rates_approx': [2.0, 0.2],  # Approximate expected rates
        'expected_weights_approx': [0.5, 0.5],  # Approximate expected weights
        'n_samples': n_samples
    }


@pytest.fixture
def mock_contact_file(tmp_path, synthetic_timeseries):
    """
    Create a mock contact pickle file for testing ParallelGibbs.
    
    Parameters
    ----------
    tmp_path : Path
        Pytest temporary directory fixture
    synthetic_timeseries : dict
        Synthetic timeseries data from fixture
        
    Returns
    -------
    str
        Path to the created contact file
    """
    from basicrta.tests.utils import make_Universe
    
    # Create simple test Universe for ag1 and ag2 with resnames and resids topology attributes
    residue_names = ['TRP', 'VAL', 'ALA', 'GLY', 'PHE', 'LEU', 'SER', 'THR', 'ASP', 'GLU']
    target_resids = [313, 314, 315, 316, 317, 318, 319, 320, 321, 322]
    
    # Create universe with topology attributes and values in one call
    ag1_universe = make_Universe(
        extras={
            'resnames': residue_names[:10],  # 10 residues
            'resids': target_resids[:10]
        }, 
        size=(50, 10, 1)  # 50 atoms, 10 residues
    )
    ag2_universe = make_Universe(size=(100, 20, 1))  # 100 atoms, 20 residues (no special attributes needed)
    
    # Create AtomGroups
    ag1 = ag1_universe.atoms
    ag2 = ag2_universe.atoms
    
    # Create contact data structure
    times = synthetic_timeseries['times']
    n_contacts = len(times)
    
    # Contact format: [protein_resid, lipid_resid, frame, residence_time, contact_number]
    contacts = np.zeros((n_contacts, 5))
    contacts[:, 0] = 313  # All contacts with residue 313 (TRP)
    contacts[:, 1] = np.arange(n_contacts)  # Different lipid residues
    contacts[:, 2] = np.arange(n_contacts)  # Sequential frames
    contacts[:, 3] = times  # Residence times
    contacts[:, 4] = np.arange(n_contacts)  # Contact numbers
    
    # Create metadata
    metadata = {
        'ag1': ag1,
        'ag2': ag2,
        'ts': 0.1,  # timestep
        'top': 'test.pdb',
        'traj': 'test.xtc'
    }
    
    # Create numpy array with metadata
    contacts_dtype = np.dtype(contacts.dtype, metadata=metadata)
    contacts_array = contacts.astype(contacts_dtype)
    
    # Save to pickle file
    contact_file = tmp_path / "contacts_7.0.pkl"
    with open(contact_file, 'wb') as f:
        pickle.dump(contacts_array, f)
        
    return str(contact_file)


class TestGibbsSampler:
    """Comprehensive tests for both Gibbs sampler classes."""

    @pytest.mark.parametrize("init_kwargs", [
        {
            'times': None,  # Will be set from fixture
            'residue': 'W313',
            'ncomp': 2,
            'niter': 1000,
            'burnin': 5,
            'cutoff': 7.0,
            'g': 100
        },
    ])
    def test_gibbs_run_method(self, tmp_path, synthetic_timeseries, init_kwargs):
        """Test the run() method for Gibbs class with synthetic data."""
        
        # Set up the times from fixture
        if 'times' in init_kwargs:
            init_kwargs['times'] = synthetic_timeseries['times']
            
        # Create output directory structure
        output_dir = tmp_path / f"basicrta-{init_kwargs['cutoff']}" / init_kwargs['residue']
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Change to tmp_path to avoid creating output in the repo
        with work_in(tmp_path):
            # Initialize Gibbs sampler
            gibbs = Gibbs(**init_kwargs)
            
            # Test that initialization worked correctly
            assert gibbs.times is not None, "Times should be set"
            assert gibbs.ncomp == 2, "Should have 2 components"
            assert gibbs.niter == 1000, "Should have 1000 iterations"
            
            # Run the sampler
            gibbs.run()
            
            # Verify that the sampler ran successfully
            assert hasattr(gibbs, 'mcweights'), "Gibbs sampler should have mcweights after running"
            assert hasattr(gibbs, 'mcrates'), "Gibbs sampler should have mcrates after running"
            assert gibbs.mcweights is not None, "mcweights should not be None"
            assert gibbs.mcrates is not None, "mcrates should not be None"
            
            # Check that we have the expected number of samples
            expected_samples = (gibbs.niter + 1) // gibbs.g
            assert gibbs.mcweights.shape[0] == expected_samples, f"Expected {expected_samples} weight samples"
            assert gibbs.mcrates.shape[0] == expected_samples, f"Expected {expected_samples} rate samples"
            
            # Check dimensions match number of components
            assert gibbs.mcweights.shape[1] == gibbs.ncomp, f"Should have {gibbs.ncomp} weight components"
            assert gibbs.mcrates.shape[1] == gibbs.ncomp, f"Should have {gibbs.ncomp} rate components"
            
            # Verify weights are properly normalized (sum to ~1 for each sample)
            weight_sums = np.sum(gibbs.mcweights, axis=1)
            np.testing.assert_allclose(weight_sums, 1.0, rtol=1e-10, 
                                     err_msg="Weights should sum to 1 for each sample")
            
            # Verify rates are positive
            assert np.all(gibbs.mcrates > 0), "All rates should be positive"
            
            # Check that survival function was computed
            assert hasattr(gibbs, 't'), "Should have time points for survival function"
            assert hasattr(gibbs, 's'), "Should have survival function values"
            assert len(gibbs.t) > 0, "Time points should not be empty"
            assert len(gibbs.s) > 0, "Survival function should not be empty"
            
            # Check that indicator variables were stored
            assert hasattr(gibbs, 'indicator'), "Should have indicator variables"
            assert gibbs.indicator is not None, "Indicator should not be None"
            assert gibbs.indicator.shape[0] == expected_samples, "Indicator should match sample count"
            assert gibbs.indicator.shape[1] == len(gibbs.times), "Indicator should match data count"


    def test_parallel_gibbs_initialization_and_contact_loading(self, tmp_path, mock_contact_file, synthetic_timeseries):
        """Test ParallelGibbs initialization and contact file loading."""
        
        # Change to tmp_path for output
        with work_in(tmp_path):
            # Initialize ParallelGibbs
            parallel_gibbs = ParallelGibbs(
                contacts=mock_contact_file,
                nproc=2,   # 2 cores are always available ...
                ncomp=4,
                niter=1000
            )
            
            # Test initialization
            assert parallel_gibbs.contacts == mock_contact_file
            assert parallel_gibbs.nproc == 2
            assert parallel_gibbs.ncomp == 4
            assert parallel_gibbs.niter == 1000
            assert parallel_gibbs.cutoff == 7.0  # Extracted from filename
            
            # Test that the contact file can be loaded and processed
            with open(mock_contact_file, 'rb') as f:
                contacts = pickle.load(f)
            
            # Verify the contact file structure
            assert contacts.shape[1] == 5, "Contact data should have 5 columns"
            assert len(contacts) == len(synthetic_timeseries['times']), "Should have correct number of contacts"
            
            # Test contact data processing logic (without running full Gibbs)
            metadata = contacts.dtype.metadata
            assert 'ag1' in metadata, "Metadata should contain ag1"
            assert 'ag2' in metadata, "Metadata should contain ag2"
            assert 'ts' in metadata, "Metadata should contain timestep"
            
            protids = np.unique(contacts[:, 0])
            assert 313 in protids, "Should have protein residue 313"
            
            # Test that residence times are extracted correctly
            residue_313_times = contacts[contacts[:, 0] == 313][:, 3]
            assert len(residue_313_times) > 0, "Should have residence times for residue 313"
            assert np.allclose(residue_313_times, synthetic_timeseries['times']), "Residence times should match synthetic data"


    def test_gibbs_initialization_parameters(self, synthetic_timeseries):
        """Test that Gibbs class initializes with correct parameters."""
        times = synthetic_timeseries['times']
        
        gibbs = Gibbs(
            times=times,
            residue='TEST123',
            loc=0,
            ncomp=3,
            niter=5000,
            cutoff=8.5,
            g=25,
            burnin=500,
            gskip=5
        )
        
        # Test all initialization parameters
        assert np.array_equal(gibbs.times, times), "Times should be stored correctly"
        assert gibbs.residue == 'TEST123', "Residue should be stored correctly"
        assert gibbs.loc == 0, "Location should be stored correctly"
        assert gibbs.ncomp == 3, "Number of components should be stored correctly"
        assert gibbs.niter == 5000, "Number of iterations should be stored correctly"
        assert gibbs.cutoff == 8.5, "Cutoff should be stored correctly"
        assert gibbs.g == 25, "Gibbs skip should be stored correctly"
        assert gibbs.burnin == 500, "Burnin should be stored correctly"
        assert gibbs.gskip == 5, "Gibbs skip should be stored correctly"
        
        # Test that timestep is computed correctly
        assert gibbs.ts is not None, "Timestep should be computed"
        assert gibbs.ts > 0, "Timestep should be positive"


    def test_gibbs_prepare_method(self, synthetic_timeseries):
        """Test the _prepare() method of Gibbs class."""
        times = synthetic_timeseries['times']
        
        gibbs = Gibbs(
            times=times,
            residue='W313',
            ncomp=2,
            niter=1000,
            g=100
        )
        
        # Call _prepare method
        gibbs._prepare()
        
        # Check that survival function was computed
        assert hasattr(gibbs, 't'), "Should have time points"
        assert hasattr(gibbs, 's'), "Should have survival function"
        assert len(gibbs.t) > 0, "Time points should not be empty"
        assert len(gibbs.s) > 0, "Survival function should not be empty"
        assert np.all(gibbs.s >= 0), "Survival function should be non-negative"
        assert np.all(gibbs.s <= 1), "Survival function should be <= 1"
        
        # Check that arrays were initialized with correct shapes
        expected_samples = (gibbs.niter + 1) // gibbs.g
        assert gibbs.indicator.shape == (expected_samples, len(times)), "Indicator shape should be correct"
        assert gibbs.mcweights.shape == (expected_samples, gibbs.ncomp), "mcweights shape should be correct"
        assert gibbs.mcrates.shape == (expected_samples, gibbs.ncomp), "mcrates shape should be correct"
        
        # Check that hyperparameters were initialized
        assert hasattr(gibbs, 'whypers'), "Should have weight hyperparameters"
        assert hasattr(gibbs, 'rhypers'), "Should have rate hyperparameters"
        assert gibbs.whypers.shape == (gibbs.ncomp,), "Weight hyperparameters should have correct shape"
        assert gibbs.rhypers.shape == (gibbs.ncomp, 2), "Rate hyperparameters should have correct shape"


    def test_parallel_gibbs_initialization(self, mock_contact_file):
        """Test ParallelGibbs initialization parameters."""
        parallel_gibbs = ParallelGibbs(
            contacts=mock_contact_file,
            nproc=4,
            ncomp=5,
            niter=50000
        )
        
        assert parallel_gibbs.contacts == mock_contact_file, "Contacts file should be stored"
        assert parallel_gibbs.nproc == 4, "Number of processes should be stored"
        assert parallel_gibbs.ncomp == 5, "Number of components should be stored"
        assert parallel_gibbs.niter == 50000, "Number of iterations should be stored"
        assert parallel_gibbs.cutoff == 7.0, "Cutoff should be extracted from filename"


    def test_parallel_gibbs_run_method(self, tmp_path, mock_contact_file, synthetic_timeseries):
        """Test the run() method for ParallelGibbs class with real multiprocessing."""
        
        # Change to tmp_path for output
        with work_in(tmp_path):
            # Initialize ParallelGibbs with smaller parameters for faster testing
            parallel_gibbs = ParallelGibbs(
                contacts=mock_contact_file,
                nproc=2,  # Use 2 processes for testing
                ncomp=2,  # Use fewer components for speed
                niter=1000  # Smaller iteration count for testing
            )
            
            # Test initialization
            assert parallel_gibbs.contacts == mock_contact_file
            assert parallel_gibbs.nproc == 2
            assert parallel_gibbs.ncomp == 2
            assert parallel_gibbs.niter == 1000
            assert parallel_gibbs.cutoff == 7.0
            
            # Run ParallelGibbs on residue 313 (which exists in our mock contact file)
            parallel_gibbs.run(run_resids=[313])
            
            # Verify that the expected output directory structure was created
            expected_residue_dir = tmp_path / f"basicrta-{parallel_gibbs.cutoff}" / "W313"
            assert expected_residue_dir.exists(), f"Residue directory should be created: {expected_residue_dir}"
            
            # Verify that the Gibbs sampler output file was created
            expected_gibbs_file = expected_residue_dir / f"gibbs_{parallel_gibbs.niter}.pkl"
            assert expected_gibbs_file.exists(), f"Gibbs output file should be created: {expected_gibbs_file}"
            
            # Load and verify the Gibbs sampler results
            from basicrta.gibbs import Gibbs
            gibbs_result = Gibbs.load(str(expected_gibbs_file))
            
            # Verify the loaded Gibbs sampler has the expected properties
            assert gibbs_result.residue == 'W313', "Residue should match"
            assert gibbs_result.ncomp == 2, "Number of components should match"
            assert gibbs_result.niter == 1000, "Number of iterations should match"
            assert gibbs_result.cutoff == 7.0, "Cutoff should match"
            
            # Verify that the Gibbs sampler ran successfully
            assert hasattr(gibbs_result, 'mcweights'), "Should have mcweights"
            assert hasattr(gibbs_result, 'mcrates'), "Should have mcrates"
            assert np.all(gibbs_result.mcweights >= 0), "mcweights should be non-negative"
            assert np.all(gibbs_result.mcrates >= 0), "mcrates should be non-negative"
            
            # Check that we have the expected number of samples
            expected_samples = (gibbs_result.niter + 1) // gibbs_result.g
            assert gibbs_result.mcweights.shape[0] == expected_samples, "Should have correct number of weight samples"
            assert gibbs_result.mcrates.shape[0] == expected_samples, "Should have correct number of rate samples"
            
            # Verify that the residence times used match our synthetic data
            assert len(gibbs_result.times) == len(synthetic_timeseries['times']), "Should use all synthetic residence times"
            assert np.allclose(gibbs_result.times, synthetic_timeseries['times']), "Residence times should match synthetic data"


    def test_gibbs_sampler_old_style_input(self, tmp_path, synthetic_timeseries):
        """Test Gibbs sampler with old-style (non-combined) input data."""
        
        # Use synthetic timeseries from fixture
        test_times = synthetic_timeseries['times']
        
        # Create temporary directory for output
        output_dir = tmp_path / "basicrta-7.0" / "test_residue"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Change to tmp_path to avoid creating output in the repo
        with work_in(tmp_path):
            # Run Gibbs sampler with old-style input
            gibbs = Gibbs(
                times=test_times,
                residue='test_residue',
                ncomp=2,  # Use 2 components for stability
                niter=1000,  # 1000 steps as requested
                burnin=5,   # 5 burnin steps as requested
                cutoff=7.0,  # Set cutoff for directory creation
                g=100  # Store samples every 100 iterations
            )
            
            # Run the sampler
            gibbs.run()
            
            # Verify that the sampler ran successfully
            assert hasattr(gibbs, 'mcweights'), "Gibbs sampler should have mcweights after running"
            assert hasattr(gibbs, 'mcrates'), "Gibbs sampler should have mcrates after running"
            assert gibbs.mcweights is not None, "mcweights should not be None"
            assert gibbs.mcrates is not None, "mcrates should not be None"
            
            # Check that we have the expected number of samples
            # The array is sized as (niter + 1) // g in _prepare()
            expected_samples = (1000 + 1) // gibbs.g
            assert gibbs.mcweights.shape[0] == expected_samples, f"Expected {expected_samples} weight samples, got {gibbs.mcweights.shape[0]}"
            assert gibbs.mcrates.shape[0] == expected_samples, f"Expected {expected_samples} rate samples, got {gibbs.mcrates.shape[0]}"
            
            # Check dimensions match number of components
            assert gibbs.mcweights.shape[1] == 2, "Should have 2 weight components"
            assert gibbs.mcrates.shape[1] == 2, "Should have 2 rate components"
            
            # Verify weights are properly normalized (sum to ~1 for each sample)
            weight_sums = np.sum(gibbs.mcweights, axis=1)
            np.testing.assert_allclose(weight_sums, 1.0, rtol=1e-10, 
                                     err_msg="Weights should sum to 1 for each sample")
            
            # Verify rates are positive
            assert np.all(gibbs.mcrates > 0), "All rates should be positive"