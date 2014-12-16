####################################################################################################
# 
# Copyright [2014] EMBL - European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in
# compliance with the License.  You may obtain a copy of
# the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software distributed under the
# License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
# either express or implied. See the License for the specific language governing permissions
# and limitations under the License.
# 
####################################################################################################

"""Utilities for standardise package"""

########################################################################

# Error types...

errors = {
    "not_built":       "RDKit could not build mol",
    "no_non_salt":     "No non-salt/solvate components",
    "multi_component": "Multiple non-salt/solvate components",
    "sanity_check":    "Molecule failed sanity check",
    "timed_out":       "Time taken shows problem with moleule"
}

class StandardiseException(Exception):

    def __init__(self, name):

        self.name = name
        self.message = errors[name]
        self.args = (errors[name], )

# StandardiseException

######

# Sanity-check for molecules
# 
# See e.g. PubChem CIDs 128221 or 20643358 for examples of things that fail

from rdkit import Chem

def sanity_check(mol):

    try:

        Chem.SanitizeMol(mol)

        Chem.MolToSmiles(mol, isomericSmiles=True)

    except ValueError as err:

        raise StandardiseException("sanity_check")

# sanity_check

######

# Time-out long running operation, as this usually indicates a problem with the molecules
# 
# http://stackoverflow.com/questions/2281850/timeout-function-if-it-takes-too-long-to-finish

from functools import wraps
import errno
import os
import signal

def timeout(seconds=2):
    
    def decorator(func):
        
        def _handle_timeout(signum, frame):
            
            raise StandardiseException("timed_out")

        def wrapper(*args, **kwargs):
            
            signal.signal(signal.SIGALRM, _handle_timeout)
            
            signal.alarm(seconds)
            
            try:

                result = func(*args, **kwargs)

            finally:

                signal.alarm(0)
                
            return result

        return wraps(func)(wrapper)

    return decorator

# timeout

########################################################################
# End
########################################################################
