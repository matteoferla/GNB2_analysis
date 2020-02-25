########################################################################################################################
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
if sys.version_info[0] < 3:
    raise EnvironmentError("Hey, caveman, use Python 3.")

__doc__ = \
    """
    PDB minimization for different GNB2 models.
    
    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2019 A.D."
__license__ = "MIT"
__version__ = "1"
__citation__ = "TBA"

########################################################################################################################


import pyrosetta
pyrosetta.init()

def relax(pose, cycles:int=2):
    scorefxn = pyrosetta.get_fa_scorefxn()
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
    relax.apply(pose)
    print('Done')

################### alpha beta gamma ######################################

native_alpha = pyrosetta.rosetta.core.pose.Pose()
pyrosetta.rosetta.core.import_pose.pose_from_file(native_alpha, 'GNB2_alpha.pdb')
relax(native_alpha, 15)#equivalent to -relax:thorough
native_alpha.dump_pdb('GNB2_alpha.r.pdb')

################### GRK2 ######################################

native_alt = pyrosetta.rosetta.core.pose.Pose()
pyrosetta.rosetta.core.import_pose.pose_from_file(native_alt, 'GNB2_GRK.pdb')
relax(native_alt, 15)
native_alt.dump_pdb('GNB2_alt.r.pdb')


################### WITHOUT alpha ######################################
import pymol2

with pymol2.PyMOL() as pymol:
    pymol.cmd.load('GNB2_alpha.pdb')
    pymol.cmd.remove('chain A')
    pymol.cmd.save('GNB2_alone.pdb')

native_alone = pyrosetta.rosetta.core.pose.Pose()
pyrosetta.rosetta.core.import_pose.pose_from_file(native_alone, 'GNB2_alone.pdb')
relax(native_alone, 15)
native_alone.dump_pdb('GNB2_alone.r.pdb')

################### phosphorylated ######################################

from Bio.SeqUtils import seq3
from michelanglo_protein import ProteinAnalyser, global_settings
global_settings.startup(data_folder='/home/matteo/Coding/Michelanglo/protein-data')
p = ProteinAnalyser(taxid=9606, uniprot='P62879').load()

native_phospho = pyrosetta.rosetta.core.pose.Pose()
pyrosetta.rosetta.core.import_pose.pose_from_file(native_phospho, 'GNB2_alone.pdb')

MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
# KinaseMover = pyrosetta.rosetta.protocols.enzymatic_movers.KinaseMover
add_variant_type_to_residue = pyrosetta.rosetta.core.pose.add_variant_type_to_residue
pose2pdb = native_phospho.pdb_info().pdb2pose
for record in p.features['PSP_modified_residues']:
    change = record['from_residue'] + '-' + record['ptm']
    if record['ptm'] == 'ub':
        continue
    elif record['ptm'] == 'p':
        patch = 'phosphorylated'
    elif record['ptm'] == 'ac':
        patch = 'acetylated'
    elif record['ptm'] == 'm1':
        patch = 'monomethylated'
    elif record['ptm'] == 'm2':
        patch = 'dimethylated'
    elif record['ptm'] == 'm3':
        patch = 'trimethylated'
    else:
        raise ValueError
    new_res = f'{seq3(p["from_residue"])}:{patch}'
    r = pose2pdb(res=int(record['residue_index']), chain='B')
    MutateResidue(target=r, new_res=new_res).apply(native_phospho)

relax(native_phospho, 15)

native_phospho.dump_pdb('GNB2_phospho.r.pdb')