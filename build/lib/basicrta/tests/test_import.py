"""
Unit and regression test for the basicrta package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import basicrta


def test_basicrta_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "basicrta" in sys.modules
