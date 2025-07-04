"""
Tests for the Gibbs sampler with old-style input.
"""

import pytest
import numpy as np
import tempfile
import os
from basicrta.gibbs import Gibbs


class TestGibbsSampler:
    """Tests for Gibbs sampler with traditional input."""

    def test_gibbs_sampler_old_style_input(self, tmp_path):
        """Test Gibbs sampler with old-style (non-combined) input data."""
        
        # Create simple synthetic test data for more reliable testing
        rng = np.random.default_rng(seed=42)
        n_samples = 100
        # Create a simple bimodal distribution
        times_short = rng.exponential(0.5, n_samples // 2)  # Fast component
        times_long = rng.exponential(5.0, n_samples // 2)   # Slow component
        test_times = np.concatenate([times_short, times_long])
        rng.shuffle(test_times)  # Mix them up
        
        # Create temporary directory for output
        output_dir = tmp_path / "basicrta-7.0" / "test_residue"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Change to tmp_path to avoid creating output in the repo
        original_cwd = os.getcwd()
        os.chdir(tmp_path)
        
        try:
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
            
        finally:
            # Restore original working directory
            os.chdir(original_cwd)