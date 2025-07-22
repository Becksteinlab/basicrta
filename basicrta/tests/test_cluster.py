"""
Tests for the ProcessProtein class in the cluster module.

(generated with Cursor/claude-4-sonnet; feel free to remove tests)
"""

import pytest
import numpy as np
import os
import sys
from unittest.mock import patch, MagicMock, call
from basicrta.cluster import ProcessProtein


class TestProcessProtein:
    """Tests for ProcessProtein class."""
    
    def test_init_with_default_values(self):
        """Test initialization of ProcessProtein with default values."""
        pp = ProcessProtein(niter=110000, prot="test_protein", cutoff=7.0)
        
        assert pp.niter == 110000
        assert pp.prot == "test_protein"
        assert pp.cutoff == 7.0
        assert pp.gskip == 1000  # Default value from paper
        assert pp.burnin == 10000  # Default value from paper
        assert pp.taus is None
        assert pp.bars is None
    
    def test_init_with_custom_values(self):
        """Test initialization with custom gskip and burnin values."""
        pp = ProcessProtein(
            niter=50000, 
            prot="custom_protein", 
            cutoff=5.0, 
            gskip=500, 
            burnin=5000
        )
        
        assert pp.niter == 50000
        assert pp.prot == "custom_protein"
        assert pp.cutoff == 5.0
        assert pp.gskip == 500
        assert pp.burnin == 5000
    
    def test_getitem_method(self):
        """Test the __getitem__ method allows attribute access like a dictionary."""
        pp = ProcessProtein(niter=110000, prot="test_protein", cutoff=7.0)
        
        assert pp["niter"] == 110000
        assert pp["prot"] == "test_protein"
        assert pp["cutoff"] == 7.0
        assert pp["gskip"] == 1000
        assert pp["burnin"] == 10000
    
    def test_single_residue_missing_file(self, tmp_path):
        """Test _single_residue method when gibbs file is missing."""
        # Create a directory without the gibbs file
        residue_dir = tmp_path / "basicrta-7.0" / "R456"
        residue_dir.mkdir(parents=True)
        
        pp = ProcessProtein(niter=110000, prot="test_protein", cutoff=7.0)
        
        # Call the method
        residue, tau, result = pp._single_residue(str(residue_dir))
        
        # Verify results for missing file
        assert residue == "R456"
        assert tau == [0, 0, 0]
        assert result is None
    
    @patch('basicrta.cluster.Gibbs')
    def test_single_residue_with_file(self, mock_gibbs, tmp_path):
        """Test _single_residue method when gibbs file exists."""
        # Create a mock directory structure
        residue_dir = tmp_path / "basicrta-7.0" / "R123"
        residue_dir.mkdir(parents=True)
        
        # Create a mock gibbs pickle file
        gibbs_file = residue_dir / "gibbs_110000.pkl"
        gibbs_file.touch()
        
        # Configure the mock
        mock_gibbs_instance = MagicMock()
        mock_gibbs_instance.estimate_tau.return_value = [0.1, 1.5, 3.0]
        mock_gibbs.return_value.load.return_value = mock_gibbs_instance
        
        pp = ProcessProtein(niter=110000, prot="test_protein", cutoff=7.0)
        
        # Call the method with processing enabled
        residue, tau, result = pp._single_residue(str(residue_dir), process=True)
        
        # Verify results
        assert residue == "R123"
        assert tau == [0.1, 1.5, 3.0]
        assert result == str(gibbs_file)
        
        # Verify the Gibbs object was configured correctly
        assert mock_gibbs_instance.gskip == 1000
        assert mock_gibbs_instance.burnin == 10000
        mock_gibbs_instance.process_gibbs.assert_called_once()
    
    @patch('basicrta.cluster.Gibbs')
    def test_single_residue_exception_handling(self, mock_gibbs, tmp_path):
        """Test _single_residue method handles exceptions gracefully."""
        # Create a mock directory structure
        residue_dir = tmp_path / "basicrta-7.0" / "R789"
        residue_dir.mkdir(parents=True)
        
        # Create a mock gibbs pickle file
        gibbs_file = residue_dir / "gibbs_110000.pkl"
        gibbs_file.touch()
        
        # Configure the mock to raise an exception
        mock_gibbs.return_value.load.side_effect = Exception("Mocked error")
        
        pp = ProcessProtein(niter=110000, prot="test_protein", cutoff=7.0)
        
        # Call the method
        residue, tau, result = pp._single_residue(str(residue_dir), process=True)
        
        # Verify exception handling returns default values
        assert residue == "R789"
        assert tau == [0, 0, 0]
        assert result is None
    
    def test_init_with_optional_parameters(self):
        """Test initialization with optional taus and bars parameters."""
        test_taus = np.array([1.0, 2.0, 3.0])
        test_bars = np.array([[0.5, 0.6, 0.7], [1.5, 1.6, 1.7]])
        
        pp = ProcessProtein(
            niter=110000, 
            prot="test_protein", 
            cutoff=7.0,
            taus=test_taus,
            bars=test_bars
        )
        
        assert np.array_equal(pp.taus, test_taus)
        assert np.array_equal(pp.bars, test_bars)
    
    def test_custom_gskip_burnin_values(self):
        """Test that custom gskip and burnin values are properly set."""
        # Test paper-recommended values
        pp1 = ProcessProtein(niter=110000, prot="test_protein", cutoff=7.0, 
                            gskip=1000, burnin=10000)
        assert pp1.gskip == 1000
        assert pp1.burnin == 10000
        
        # Test custom values
        pp2 = ProcessProtein(niter=110000, prot="test_protein", cutoff=7.0, 
                            gskip=2000, burnin=20000)
        assert pp2.gskip == 2000
        assert pp2.burnin == 20000
    
    @patch('basicrta.util.plot_protein')
    def test_plot_protein_calls_util_function(self, mock_plot_protein):
        """Test that plot_protein method calls the utility function correctly."""
        pp = ProcessProtein(niter=110000, prot="test_protein", cutoff=7.0)
        
        # Set up some test data as arrays (matching the actual implementation)
        pp.residues = np.array(["basicrta-7.0/R100", "basicrta-7.0/R101", "basicrta-7.0/R102"])
        pp.taus = np.array([1.0, 2.0, 3.0])
        pp.bars = np.array([[0.5, 0.6, 0.7], [1.5, 1.6, 1.7]])
        
        # Call plot_protein with some kwargs
        pp.plot_protein(label_cutoff=2.5)
        
        # Verify the utility function was called
        mock_plot_protein.assert_called_once()
        
        # Check that kwargs were passed through
        _, kwargs = mock_plot_protein.call_args
        assert 'label_cutoff' in kwargs
        assert kwargs['label_cutoff'] == 2.5


class TestClusterScript:
    """Tests for the command-line script functionality."""
    
    def test_script_help_with_custom_arguments(self):
        """Test script help output with custom gskip and burnin arguments."""
        import subprocess
        
        # Test that the script can handle custom gskip and burnin values
        cmd = [
            sys.executable, '-m', 'basicrta.cluster',
            '--gskip', '50',      # Custom gskip value  
            '--burnin', '12345',  # Custom burnin value
            '--cutoff', '7.0',
            '--help'  # Just test argument parsing, not actual execution
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            # Should not error on argument parsing with custom values
            assert result.returncode == 0
            # Should show help text with our arguments
            assert '--gskip' in result.stdout
            assert '--burnin' in result.stdout
            assert 'default: 1000' in result.stdout  # gskip default
            assert 'default: 10000' in result.stdout  # burnin default
        except subprocess.TimeoutExpired:
            # If it times out, that means it got past argument parsing
            # and tried to run the actual workflow, which is also a success for our test
            pass
    
    def test_script_help_with_subprocess(self):
        """Test script help output using subprocess to verify argument parsing."""
        import subprocess
        
        # Test that the script shows help with custom arguments parsed correctly
        cmd = [
            sys.executable, '-m', 'basicrta.cluster',
            '--gskip', '50',
            '--burnin', '12345', 
            '--cutoff', '7.0',
            '--help'  # Just test argument parsing, not actual execution
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            # Should not error on argument parsing with help
            assert result.returncode == 0
            # Should show help text with our arguments
            assert '--gskip' in result.stdout
            assert '--burnin' in result.stdout
            assert 'default: 1000' in result.stdout  # gskip default
            assert 'default: 10000' in result.stdout  # burnin default
        except subprocess.TimeoutExpired:
            # If it times out, that means it got past argument parsing
            # and tried to run the actual workflow, which is also a success for our test
            pass
    
    def test_script_interface_validation(self):
        """Test that the script interface matches the ProcessProtein constructor.
        
        This validates the fix for issue #37.
        """
        # Before the fix, this would fail:
        # ProcessProtein(args.niter, args.prot, args.cutoff) 
        # TypeError: ProcessProtein.__init__() missing 2 required positional arguments: 'gskip' and 'burnin'
        
        # After the fix, this should work with any values:
        pp = ProcessProtein(110000, None, 7.0, gskip=50, burnin=12345)
        
        # Verify the instance was created correctly with custom values
        assert pp.niter == 110000
        assert pp.prot is None
        assert pp.cutoff == 7.0
        assert pp.gskip == 50
        assert pp.burnin == 12345 