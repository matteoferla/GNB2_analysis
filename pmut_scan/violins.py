########################################################################################################################
# !/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

if sys.version_info[0] < 3:
    raise EnvironmentError("Hey, caveman, use Python 3.")

__doc__ = \
    """
    The following code was actually, in a Jupyter notebook so it may not work.
    Also includes file paths.

    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2019 A.D."
__license__ = "MIT"
__version__ = "1"
__citation__ = "TBA"

########################################################################################################################

import re
import plotly.graph_objects as go
import pandas as pd

def read_data(grouping: str) -> pd.DataFrame:
    data = []
    with open('/home/matteo/Desktop/GNB2_part2/models_try-1/' + grouping + '.txt') as fh:  # unpaired
        for line in fh:
            line = line.replace('protocols.pmut_scan.PointMutScanDriver: ', '')
            if not line.strip() or '()' in line or 'mutation' in line:
                continue
            d = dict(zip(['xmutation', 'mutation', 'ddG', 'total_energy'], line.split()))
            if 'total_energy' not in d:
                continue  # len check basically
            rex = re.match('(\w)-(\w)(\d+)(\w)', d['mutation'])
            d['chain'] = rex.group(1)
            d['from_resn'] = rex.group(2)
            d['resi'] = int(rex.group(3))
            d['to_resn'] = rex.group(4)
            d['total_energy'] = float(d['total_energy'])
            d['ddG'] = float(d['ddG'])
            if d['chain'] != 'B':
                continue
            data.append(d)
    df = pd.DataFrame(data)
    return df.assign(mutation=df.mutation.apply(lambda v: v[2:]))

#######################################################################################################################

def reg_violins():
    for grouping, title in (('unpaired', 'ΔΔG landscape of Gβ2 (GNB2) bound to Gγ1 (GNG1) subunit only'),
                            ('wAlpha', 'ΔΔG landscape of Gβ2 bound to γ1 and αi (GNAI) subunits'),
                            ('wGRK2', 'ΔΔG landscape of Gβ2 bound to γ1 and β-adrenergic receptor kinase 1 (GRK2))')
                            ):
        df = read_data(grouping)
        df = df.assign(ddG_limited=df.ddG.apply(lambda v: min(abs(v), 50) * abs(v) / (v + 0.0001)))

        fig = go.Figure()
        for groupname, group in (('gnomAD', gnomad), ('Pathogenic', pathogenic), ('Clinvar_homologues', clinvar)):
            fig.add_trace(go.Violin(y=df.loc[df.mutation.isin(group)].ddG_limited,
                                    name=groupname,
                                    box_visible=True,
                                    meanline_visible=True))

        fig.add_trace(go.Violin(y=df.ddG_limited,
                                name='sequence-space',
                                box_visible=True,
                                meanline_visible=True))

        fig.update_layout(title_text=title, yaxis={'range': [-10, 51]})
        fig.write_image('violin_' + grouping + '.png', scale=3)
        fig.show()

#######################################################################################################################

def get_get_difference():
    ref = read_data('unpaired')
    refdex = dict(zip(ref.mutation, ref.ddG))

    def get_difference(row):
        if row.mutation in refdex:
            return row.ddG - refdex[row.mutation]
        else:
            return 0

    return get_difference


def diff_violins():
    get_difference = get_get_difference()
    for grouping, title in (
            ('wAlpha', 'ΔΔG landscape of Gβ2 bound to γ1 and αi (GNAI) subunits'),
            ('wGRK2', 'ΔΔG landscape of Gβ2 bound to γ1 and β-adrenergic receptor kinase 1 (GRK2))')
    ):

        df = read_data(grouping)
        df = df.assign(dddG=df.apply(get_difference, 1))
        df = df.assign(dddG_limited=df.dddG.apply(lambda v: min(abs(v), 50) * abs(v) / (v + 0.0001)))

        fig = go.Figure()
        for groupname, group in (('gnomAD', gnomad), ('Pathogenic', pathogenic), ('Clinvar_homologues', clinvar)):
            fig.add_trace(go.Violin(y=df.loc[df.mutation.isin(group)].dddG_limited,
                                    name=groupname,
                                    box_visible=True,
                                    meanline_visible=True))
        fig.add_trace(go.Violin(y=df.dddG_limited,
                                name='sequence-space',
                                box_visible=True,
                                meanline_visible=True))

        fig.update_layout(title_text=title, yaxis={'range': [-10, 51]})
        fig.write_image('violin_diff_' + grouping + '.png', scale=3)
        fig.show()

#######################################################################################################################

reg_violins()
diff_violins()
