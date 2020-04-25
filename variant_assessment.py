########################################################################################################################
# !/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

if sys.version_info[0] < 3:
    raise EnvironmentError("Hey, caveman, use Python 3.")

__doc__ = \
    """
    ddG calculation for variants. Mainly salvaged from my DogCatcher repo.
    In variant assessment old is an added interface scoring. But I do not trust the repacking.

    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2020 A.D."
__license__ = "MIT"
__version__ = "1"
__citation__ = "TBA"

########################################################################################################################

import re, csv, os
from typing import Dict
import pyrosetta

pyrosetta.init(extra_options='-mute all')



class Variant:
    """
    Copy pasted from DogCatcher.
    """

    _name3 = {'A': 'ALA',
             'C': 'CYS',
             'D': 'ASP',
             'E': 'GLU',
             'F': 'PHE',
             'G': 'GLY',
             'H': 'HIS',
             'I': 'ILE',
             'L': 'LEU',
             'K': 'LYS',
             'M': 'MET',
             'N': 'ASN',
             'P': 'PRO',
             'Q': 'GLN',
             'R': 'ARG',
             'S': 'SER',
             'T': 'THR',
             'V': 'VAL',
             'W': 'TRP',
             'Y': 'TYR'}

    def __init__(self, filename:str):
        self.pose = self.load_pose_from_file(filename)

    def load_pose_from_file(self, filename: str) -> pyrosetta.Pose:
        """
        Loads a pose from filename with the params in the params_folder

        :param filename:
        :return:
        """
        pose = pyrosetta.Pose()
        # params_paths = pyrosetta.rosetta.utility.vector1_string()
        # params_paths.extend([os.path.join(self.params_folder, file) for file in os.listdir(self.params_folder)
        #                      if os.path.splitext(file)[1] == '.params'])
        # pyrosetta.generate_nonstandard_residue_set(pose, params_paths)
        pyrosetta.rosetta.core.import_pose.pose_from_file(pose, filename)
        return pose

    def relax_around_mover(self,
                           pose: pyrosetta.Pose,
                           resi: int, chain: str,
                           scorefxn=None, cycles=5, distance=5, cartesian=False) -> None:
        """
        Relaxes pose ``distance`` around resi:chain.

        :param resi: PDB residue number.
        :param chain:
        :param pose:
        :param scorefxn:
        :param cycles: of relax (3 quick, 15 thorough)
        :param distance:
        :param cartesian:
        :return:
        """
        if scorefxn is None:
            scorefxn = pyrosetta.get_fa_scorefxn()
            #self._cst_score(scorefxn)
        movemap = pyrosetta.MoveMap()
        ####
        resi_sele = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
        resi_sele.set_index(pose.pdb_info().pdb2pose(chain=chain, res=resi))
        NeighborhoodResidueSelector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector
        neigh_sele = NeighborhoodResidueSelector(resi_sele, distance=distance, include_focus_in_subset=True)
        n = neigh_sele.apply(pose)
        movemap.set_bb(allow_bb=n)
        movemap.set_chi(allow_chi=n)
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
        relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
        relax.set_movemap(movemap)
        relax.cartesian(cartesian)
        relax.apply(pose)

    def make_mutant(self, pose: pyrosetta.Pose, mutation:str, chain='A') -> pyrosetta.Pose:
        """
        Make a point mutant (``A23D``).

        :param pose: pose
        :param mutation:
        :param chain:
        :return:
        """
        mutant = pose.clone()
        pose2pdb = pose.pdb_info().pdb2pose
        rex = re.match('(\w)(\d+)(\w)', mutation)
        r = pose2pdb(res=int(rex.group(2)), chain=chain)
        rn = pose.residue(r).name1()
        assert rn == rex.group(1), f'residue {r}(pose)/{rex.group(2)}(pdb) is a {rn}, not a {rex.group()}'
        MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
        MutateResidue(target=r, new_res=self._name3[rex.group(3)]).apply(mutant)
        self.relax_around_mover(mutant, int(rex.group(2)), chain, distance=12, cycles=15)
        return mutant


if __name__ == '__main__':
    with open('data/GNB2_analysis2x.csv', 'w') as w:
        out = csv.DictWriter(w, fieldnames=['model',
                                            'mutation',
                                            'ddG',
                                            'native_score',
                                            'mutant_score'])
        out.writeheader()
        scorefxn = pyrosetta.get_fa_scorefxn()
        for filename in os.listdir('templates'):
            if '.pdb' not in filename:
                continue
            if filename != 'GNB2_wAlpha.pdb':
                continue
            print(filename)
            modelname = re.search(r"GNB2_(.*)\.pdb", filename).group(1)
            model = Variant(os.path.join('templates',filename))
            for mutation in ('R52L', 'A73T', 'G77R', 'G77E', 'K89T', 'K89E', 'S147L'):
                print(mutation)
                variant = model.make_mutant(model.pose, mutation=mutation, chain='B')
                variant.dump_pdb(f'variants/{modelname}.{mutation}.pdb')
                n = scorefxn(model.pose)
                m = scorefxn(variant)
                out.writerow({'model': modelname,
                              'mutation': mutation,
                              'ddG': m - n,
                              'native_score': n,
                              'mutant_score': m})




