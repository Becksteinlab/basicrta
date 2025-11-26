"""
Tests for combining contact timeseries from multiple repeat runs.
"""
import basicrta
import os
import pytest
import numpy as np
import pickle
import subprocess
import basicrta.cli
import importlib
import argparse
from basicrta.contacts import CombineContacts


class TestCLI:
    """Test class for cli.py functionality."""
    modules = basicrta.cli.commands
        
    def test_cli_modules(self):
        """Test cli as module"""
        for module in self.modules:
            #help successfully printed
            help_ret = subprocess.run(['python', '-m', f'basicrta.{module}', 
                                     '--help'])
            assert help_ret.returncode == 0
            
            # error if required arguments not given
            nohelp = subprocess.run(['python', '-m', f'basicrta.{module}'])
            assert nohelp.returncode == 2

    def test_cli_entrypoint(self):
        # print general help if no command given
        assert subprocess.run('basicrta').returncode == 0

        for module in self.modules:
            #help successfully printed
            help_ret = subprocess.run(['basicrta', f'{module}', '--help'])
            assert help_ret.returncode == 0

            # help is printed if no arguments given
            nohelp = subprocess.run(['basicrta', f'{module}'])
            assert nohelp.returncode == 0

    def test_get_module_parsers(self):
        for module in self.modules:
            parser = importlib.import_module(f"basicrta.{module}").get_parser()
            assert type(parser) == argparse.ArgumentParser

    def test_call_main_empty(self):
        for module in self.modules:
            with pytest.raises(SystemExit):
                importlib.import_module(f"basicrta.{module}").main()

    def test_cli_script_call(self):
            #help successfully printed
            help_ret = subprocess.run(['python', '-m', 'basicrta.cli', 
                                       '--help'])
            assert help_ret.returncode == 0
            
            # error if required arguments not given
            nohelp = subprocess.run(['python', '-m', 'basicrta.cli'])
            assert nohelp.returncode == 2


#    def test_call_main_args(self):
#        with mock.patch('sys.argv', ['cluster', '--cutoff', '6.9']):
#            importlib.import_module(f"basicrta.cluster").main()


