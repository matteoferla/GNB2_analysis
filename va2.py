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
    inspired from https://graylab.jhu.edu/pyrosetta/downloads/scripts/demo/D090_Ala_scan.py

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
        inspired from https://graylab.jhu.edu/pyrosetta/downloads/scripts/demo/D090_Ala_scan.py

        :param pose: changed in place
        :param chain: chain
        :return:
        """
    pdb2pose = pose.pdb_info().pdb2pose
    for i in range(1, 999):
        r = pdb2pose(res=i, chain=chain)
        if r != 0:
            shift_residue(pose, r)


def repack_interface(pose: pyrosetta.rosetta.core.pose.Pose, chain) -> None:
    scorefxn = pyrosetta.get_fa_scorefxn()
    # packer_task = pyrosetta.rosetta.core.pack.task.TaskFactory.create_packer_task(pose)
    # packer_task.restrict_to_repacking()
    # pyrosetta.rosetta.core.pack.pack_rotamers(pose, scorefxn, packer_task)
    operation = pyrosetta.rosetta.core.pack.task.operation
    allow = operation.RestrictToRepackingRLT()
    chainB_sele = pyrosetta.rosetta.core.select.residue_selector.ChainSelector('B')
    not_chainB_sele = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(chainB_sele)
    sele = pyrosetta.rosetta.core.select.residue_selector.InterGroupInterfaceByVectorSelector(chainB_sele, not_chainB_sele)
    sele.nearby_atom_cut(6)
    sele.cb_dist_cut(8)
    restrict_to_focus = operation.OperateOnResidueSubset(allow, sele, True)
    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
    tf.push_back(operation.PreventRepacking())
    tf.push_back(restrict_to_focus)
    packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn)
    packer.task_factory(tf)
    packer.apply(pose)
    print('done')


def interface_score(pose: pyrosetta.rosetta.core.pose.Pose, chain: str) -> Dict:
    # Copy
    shift_pose = pyrosetta.rosetta.core.pose.Pose()
    shift_pose.assign(pose)
    # Shift
    shift_chain(shift_pose, chain)
    repack_interface(shift_pose, chain)
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

global_settings.startup(data_folder=os.environ['MICHELANGLODATA'])

from michelanglo_protein.analyse import Mutator


def get_data(uniprot: str) -> ProteinCore:
    return ProteinCore(taxid=9606, uniprot=uniprot).load()


def assess_variant(mutation: str, groupname: str, model: str, pdbblock: str) -> Dict:
    from_res = mutation[0]
    resi = int(re.search('^\w(\d+)\w', mutation).group(1))
    to_res = mutation[-1]
    # lazy reuse of code.
    m = Mutator(pdbblock=pdbblock,
                target_resi=resi,
                target_chain='B',
                cycles=15,
                radius=5)
    analysis = m.analyse_mutation(to_res)
    m.pose.dump_pdb(f'variants/{model}.{mutation}.pdb')
    shifted = interface_score(m.pose, 'B')
    analysis = {**analysis, **shifted}
    ref = interface_score(m.native, 'B')
    shifted['apo_pose'].dump_pdb(f'variants/{model}.{mutation}.split.pdb')
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
    #this code requires Michelanglo data, solely to get the gnomads...
    #p = get_data('P62879')
    #gnomad = [re.match('^(\w\d+\w)', variant.description).group(1) for variant in p.gnomAD]
    gnomad = [] #['A21S', 'R49K', 'A73T', 'R129H', 'V133I', 'D195N', 'R197H', 'T198M', 'I208V', 'V213M', 'M263V', 'R283W', 'G302S', 'D303N', 'A309T', 'D322H', 'D323N', 'S2R', 'E12D', 'R22Q', 'T31I', 'G36R', 'G36E', 'D38E', 'I43M', 'R49K', 'D66E', 'V71L', 'L79V', 'S84T', 'N88T', 'A104T', 'Y105F', 'C114S', 'I120V', 'I123V', 'I123T', 'R129H', 'R137G', 'T143A', 'I157V', 'T159S', 'T164A', 'T173I', 'G174S', 'V178L', 'S191A', 'A193S', 'D195N', 'R197C', 'R197H', 'T198K', 'S207C', 'I208V', 'K209R', 'D212H', 'V213M', 'M217V', 'R219Q', 'I223V', 'F241V', 'G244V', 'A248V', 'T249M', 'F253L', 'D258Y', 'L262V', 'M263V', 'M263T', 'H266N', 'H266R', 'N268D', 'N268I', 'G272S', 'S275A', 'R280C', 'R280H', 'A287T', 'I296T', 'A299T', 'M300V', 'G302S', 'D303N', 'R304H', 'A305T', 'A309G', 'D312V', 'V315M', 'L318I', 'G319X', 'D322N', 'D322G', 'D323N', 'M325V', 'V327M', 'F335L', 'I338V']
    pathogenic = ('R52L',) #('A73T', 'G77R', 'G77E', 'K89T', 'K89E', 'E180K', 'S147L', 'I171T')
    table = []

    with open('data/GNB2_analysis.csv', 'a') as w:
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
            model = re.search(r"GNB2_(.*)\.pdb", filename).group(1)
            pdbblock = open(os.path.join('templates', filename)).read()
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


if __name__ == '__main__':
    main()
