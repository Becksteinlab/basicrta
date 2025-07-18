# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# Basicrta
# Copyright (c) 2024 Ricky Sexton
#
# Released under the GNU Public Licence, v3 or any higher version
#
# Please cite your use of basicrta in published work:
#

"""
Location of data files for the basicrta unit tests
====================================================

  from basicrta.tests.datafiles import *

"""

__all__ = [
    "times" # set of residence times collected from W313 of b2ar 
]

from importlib import resources
import numpy as np
import basicrta.tests.data

_data_ref = resources.files('basicrta.tests.data')

times = np.load((_data_ref / 'times.npy').as_posix())
PDB = (_data_ref / 'prot_chol.pdb').as_posix()
XTC = (_data_ref / 'prot_chol.xtc').as_posix()
# This should be the last line: clean up namespace
del resources

