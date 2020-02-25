########################################################################################################################
# !/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

if sys.version_info[0] < 3:
    raise EnvironmentError("Hey, caveman, use Python 3.")

__doc__ = \
    """
    ddG calculation for variants.

    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2019 A.D."
__license__ = "MIT"
__version__ = "1"
__citation__ = "TBA"

########################################################################################################################

import re, csv, os
from typing import Dict
import pyrosetta
pyrosetta.init()


def relax(pose, cycles: int = 2) -> None:
    """
    Relax

    :param pose: changed in place
    :param cycles:
    :return:
    """
    scorefxn = pyrosetta.get_fa_scorefxn()
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
    relax.apply(pose)
    print('Done')


def shift_residue(pose: pyrosetta.rosetta.core.pose.Pose, r: int) -> None:
    """
    copied from https://graylab.jhu.edu/pyrosetta/downloads/scripts/demo/D090_Ala_scan.py

    :param pose: changed in place
    :param r: pose numbering
    :return:
    """
    xyz = pyrosetta.rosetta.numeric.xyzVector_double_t()
    xyz.x = 500.0
    xyz.y = 0.0
    xyz.z = 0.0
    for a in range(1, pose.residue(r).natoms() + 1):
        pose.residue(r).set_xyz(a,
                                pose.residue(r).xyz(a) + xyz)


def shift_chain(pose: pyrosetta.rosetta.core.pose.Pose, chain: str) -> None:
    """
        copied from https://graylab.jhu.edu/pyrosetta/downloads/scripts/demo/D090_Ala_scan.py

        :param pose: changed in place
        :param chain: chain
        :return:
        """
    pdb2pose = pose.pdb_info().pdb2pose
    for i in range(1, 999):
        r = pdb2pose(res=i, chain=chain)
        if r != 0:
            shift_residue(pose, r)


def repack(pose: pyrosetta.rosetta.core.pose.Pose) -> None:
    scorefxn = pyrosetta.get_fa_scorefxn()
    packer_task = pyrosetta.rosetta.core.pack.task.TaskFactory.create_packer_task(pose)
    packer_task.restrict_to_repacking()
    pyrosetta.rosetta.core.pack.pack_rotamers(pose, scorefxn, packer_task)
    print('done')


def interface_score(pose: pyrosetta.rosetta.core.pose.Pose, chain: str) -> Dict:
    # Copy
    shift_pose = pyrosetta.rosetta.core.pose.Pose()
    shift_pose.assign(pose)
    # Shift
    shift_chain(shift_pose, chain)
    repack(shift_pose)
    # Score
    scorefxn = pyrosetta.get_fa_scorefxn()
    holoscore = scorefxn(pose)
    aposcore = scorefxn(shift_pose)
    #     print('holoscore', holoscore)
    #     print('aposcore', aposcore)
    #     shift_pose.dump_pdb('shifted.pdb')
    #     pose.dump_pdb('unshifted.pdb')
    #     raise ValueError
    return {'holo_score': holoscore,
            'apo_score': aposcore,
            'interface': holoscore - aposcore,
            'holo_pose': pose,
            'apo_pose': shift_pose}

########################################################################################################################

from michelanglo_protein import ProteinCore, global_settings
global_settings.startup(data_folder='/home/matteo/Coding/Michelanglo/protein-data')

from michelanglo_protein.analyse import Mutator

def get_data(uniprot:str) -> ProteinCore:
    p = ProteinCore(taxid=9606, uniprot=uniprot).load()

def assess_variant(mutation: str, groupname: str, model: str, pdbblock: str) -> Dict:
    from_res = mutation[0]
    resi = int(re.search('^\w(\d+)\w', mutation).group(1))
    to_res = mutation[-1]
    m = Mutator(pdbblock=pdbblock,
                target_resi=resi,
                target_chain='B',
                cycles=3,
                radius=4)
    analysis = m.analyse_mutation(to_res)
    m.pose.dump_pdb(f'mutants/{model}.{mutation}.pdb')
    shifted = interface_score(m.pose, 'B')
    analysis = {**analysis, **shifted}
    ref = interface_score(m.native, 'B')
    shifted['apo_pose'].dump_pdb(f'mutants/{model}.{mutation}.split.pdb')
    return {
        'model': model,
        'groupname': groupname,
        'mutation': mutation,
        'ddG': analysis['ddG'],
        'native_score': analysis['scores']['relaxed'],
        'mutant_score_fixed': analysis['scores']['mutate'],
        'mutant_score': analysis['scores']['mutarelax'],
        'interface_score': analysis['interface'],
        'interface_Δscore': analysis['interface'] - ref['interface'],
        'rmsd': analysis['rmsd']
    }

def main():
    p = get_data('P62879')
    gnomad = [re.match('^(\w\d+\w)', variant.description).group(1) for variant in p.gnomAD]
    pathogenic = ('A73T', 'G77R', 'G77E', 'K89T', 'K89E', 'E180K', 'S147L', 'I171T')
    table = []

    with open('GNB2_analysis.csv','w') as w:
        out = csv.DictWriter(w, fieldnames=['model',
                                            'groupname',
                                            'mutation',
                                            'ddG',
                                            'native_score',
                                            'mutant_score_fixed',
                                            'mutant_score',
                                            'interface_score',
                                            'interface_Δscore',
                                           'rmsd'])
        out.writeheader()
        for filename in os.listdir('templates'):
            if '.pdb' not in filename:
                continue
            print(filename)
            model = re.search(r"GNB2_(.*)\.pdb",filename).group(1)
            pdbblock = open(os.path.join('GNB2',filename)).read()
            for group, groupname in [(pathogenic, 'pathogenic'), (gnomad, 'gnomad')]:
                for mutation in group:
                    if 'X' in mutation:
                        continue
                    if 'fs' in mutation:
                        continue
                    if '*' in mutation:
                        continue
                    entry = assess_variant(mutation, groupname, model, pdbblock)
                    table.append(entry)
                    out.writerow(entry)