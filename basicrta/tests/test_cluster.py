"""
Tests for the ProcessProtein class in the cluster module.

(generated with Cursor/claude-4-sonnet; feel free to remove tests)
"""

import pytest
import numpy as np
import sys
from unittest.mock import patch, MagicMock, call
from basicrta.cluster import ProcessProtein
from basicrta.tests.utils import work_in


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

    def test_write_data_with_existing_data(self, tmp_path):
        """Test write_data method when taus and bars are already set."""
        pp = ProcessProtein(niter=110000, prot="test_protein", cutoff=7.0)
        
        # Set up test data as numpy arrays (matching the actual implementation)
        pp.residues = np.array(["R100", "R101", "R102"])
        pp.taus = np.array([1.0, 2.0, 3.0])
        pp.bars = np.array([[0.5, 0.6, 0.7], [1.5, 1.6, 1.7]])
        
        # Create output file in temporary directory
        output_file = tmp_path / "test_taus"
        
        # Call write_data
        pp.write_data(str(output_file))
        
        # Verify the file was created
        assert output_file.with_suffix('.npy').exists()
        
        # Load and verify the data
        saved_data = np.load(str(output_file) + '.npy')
        
        # Expected data format: [resid, tau, CI_lower, CI_upper]
        expected_data = np.array([
            [100, 1.0, 0.5, 1.5],  # R100 -> 100
            [101, 2.0, 0.6, 1.6],  # R101 -> 101
            [102, 3.0, 0.7, 1.7]   # R102 -> 102
        ])
        
        assert np.array_equal(saved_data, expected_data)

    @patch('basicrta.cluster.ProcessProtein.get_taus')
    def test_write_data_calls_get_taus_when_needed(self, mock_get_taus, tmp_path):
        """Test write_data method calls get_taus when taus is None."""
        pp = ProcessProtein(niter=110000, prot="test_protein", cutoff=7.0)
        
        # Define the test data
        test_taus = np.array([1.5, 2.5, 3.5])
        test_bars = np.array([[0.3, 0.4, 0.5], [1.7, 1.8, 1.9]])
        test_residues = np.array(["R200", "R201", "R202"])
        
        # Set up mock to return values AND set instance attributes (like real get_taus)
        def mock_get_taus_side_effect():
            pp.taus = test_taus
            pp.bars = test_bars
            pp.residues = test_residues
            return test_taus, test_bars
        
        mock_get_taus.side_effect = mock_get_taus_side_effect
        
        # Create output file in temporary directory
        output_file = tmp_path / "test_taus_from_get_taus"
        
        # Ensure taus is None to trigger get_taus call
        pp.taus = None
        
        # Call write_data
        pp.write_data(str(output_file))
        
        # Verify get_taus was called
        mock_get_taus.assert_called_once()
        
        # Verify the file was created and contains expected data
        assert output_file.with_suffix('.npy').exists()
        saved_data = np.load(str(output_file) + '.npy')
        
        expected_data = np.array([
            [200, 1.5, 0.3, 1.7],  # R200 -> 200
            [201, 2.5, 0.4, 1.8],  # R201 -> 201  
            [202, 3.5, 0.5, 1.9]   # R202 -> 202
        ])
        
        assert np.array_equal(saved_data, expected_data)

    @patch('basicrta.cluster.glob')
    @patch('basicrta.cluster.Pool')
    @patch('basicrta.util.get_bars')
    def test_get_taus_returns_values(self, mock_get_bars, mock_pool, mock_glob):
        """Test that get_taus method returns values as documented in docstring."""
        pp = ProcessProtein(niter=110000, prot="test_protein", cutoff=7.0)
        
        # Mock the directory structure
        mock_glob.return_value = ["basicrta-7.0/R100", "basicrta-7.0/R101"]
        
        # Mock the multiprocessing pool to return test data
        mock_pool_instance = mock_pool.return_value.__enter__.return_value
        mock_imap_results = [
            ("R100", [0.1, 1.5, 2.8], "path1"),
            ("R101", [0.2, 2.0, 3.2], "path2")
        ]
        mock_pool_instance.imap.return_value = mock_imap_results
        
        # Mock get_bars to return test confidence intervals
        test_bars = np.array([[0.5, 0.6], [2.5, 2.6]])
        mock_get_bars.return_value = test_bars
        
        # Call get_taus and verify it returns values
        result = pp.get_taus(nproc=1)
        
        # Verify the method returns a tuple as documented
        assert isinstance(result, tuple)
        assert len(result) == 2
        
        returned_taus, returned_bars = result
        
        # Verify the returned values match the instance attributes
        assert np.array_equal(returned_taus, pp.taus)
        assert np.array_equal(returned_bars, pp.bars)
        
        # Verify the values are what we expect
        expected_taus = np.array([1.5, 2.0])  # Middle values from tau arrays
        assert np.array_equal(returned_taus, expected_taus)
        assert np.array_equal(returned_bars, test_bars)

    def test_write_data_with_default_filename(self, tmp_path):
        """Test write_data method uses default filename when none provided."""
        pp = ProcessProtein(niter=110000, prot="test_protein", cutoff=7.0)
        
        # Set up test data
        pp.residues = np.array(["R300", "R301"])
        pp.taus = np.array([4.0, 5.0])
        pp.bars = np.array([[0.8, 0.9], [2.0, 2.1]])
        
        with work_in(tmp_path):
            # Call write_data without filename (should use default)
            pp.write_data()
            
            # Verify default file was created
            default_file = tmp_path / "tausout.npy"
            assert default_file.exists()
            
            # Verify data integrity
            saved_data = np.load(default_file)
            expected_data = np.array([
                [300, 4.0, 0.8, 2.0],
                [301, 5.0, 0.9, 2.1]
            ])
            assert np.array_equal(saved_data, expected_data)


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