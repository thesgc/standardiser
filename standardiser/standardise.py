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

"""Apply standardisation procedure"""

########################################################################

from __future__ import print_function, division

import logging; logger = logging.getLogger(__name__)

from rdkit import Chem

from standardiser import break_bonds, neutralise, rules, unsalt

from standardiser.utils import StandardiseException, sanity_check, timeout

########################################################################
#
# Module configuration...
#

########################################################################
#
# Module initialization...
#

########################################################################

#@timeout()
def apply(input_mol, output_rules_applied=None): 
    # Get input molecule...

    if type(input_mol) == Chem.rdchem.Mol:

        mol = input_mol

        input_type = 'mol'

    else:

        mol = Chem.MolFromMolBlock(input_mol)

        if not mol:

            mol = Chem.MolFromSmiles(input_mol)

            if not mol:

                raise StandardiseException("not_built")

            else:

                input_type = 'smi'
        else:

            input_type = 'sdf'

    sanity_check(mol)

    ######

    # Get disconnected fragments...

    non_salt_frags = []

    mol = break_bonds.apply(mol)
    inorganic_count = 0
    salt_count = 0
    for n, frag in enumerate(Chem.GetMolFrags(mol, asMols=True), 1):
        thispassrules = []
        changes = {}
        logger.debug("Starting fragment {n} '{smi}'...".format(n=n, smi=Chem.MolToSmiles(frag)))

        logger.debug("1) Check for non-organic elements...")

        if unsalt.is_nonorganic(frag): 
            inorganic_count += 1
            continue

        logger.debug("2) Attempting to neutralise (first pass)...")

        neutfrag = neutralise.apply(frag)
        neutralised = False
        if neutfrag != frag:
            neutralised = True

        logger.debug("3) Applying rules...")

        tautfrag = rules.apply(neutfrag, output_rules_applied=thispassrules)


        logger.debug("4) Attempting to neutralise (second pass)...")

        neutfrag2 = neutralise.apply(tautfrag)

        if neutfrag2 != tautfrag:
            neutralised = True


        logger.debug("5) Checking if frag is a salt/solvate...")

        if unsalt.is_salt(frag): 
            salt_count += 1
            continue

        logger.debug("...fragment kept.")


        non_salt_frags.append(frag)
        output_rules_applied.extend(thispassrules)

    if len(non_salt_frags) == 0:
        output_rules_applied.append(("no_non_salt", "ERROR"))

        # raise StandardiseException("no_non_salt")

    if len(non_salt_frags) > 1:
        output_rules_applied.append(("multi_component", "ERROR"))
        # raise StandardiseException("multi_component")


    parent = non_salt_frags[0]
    
    output_rules_applied.append(("inorganic_count", inorganic_count))
    output_rules_applied.append(("salt_count", salt_count))
    ######

    # Return parent in same format as input...

    if input_type == 'mol':

        return parent

    elif input_type == 'sdf':

        return Chem.MolToMolBlock(parent)

    else:  # input_type == 'smi'

        return Chem.MolToSmiles(parent, isomericSmiles=True)

# apply

########################################################################
# End
########################################################################
