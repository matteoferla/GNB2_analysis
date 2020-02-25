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
import pandas as pd
import plotly.graph_objects as go

this_chain = 'B'
omni = {}
for grouping, title in (('unpaired', 'ΔΔG landscape of Gβ2 (GNB2) bound to Gγ1 (GNG1) subunit only'),
                        ('wAlpha', 'ΔΔG landscape of Gβ2 bound to γ1 and αi (GNAI) subunits'),
                        ('wGRK2', 'ΔΔG landscape of Gβ2 bound to γ1 and β-adrenergic receptor kinase 1 (GRK2))')
                        ):
    # protocols.pmut_scan.PointMutScanDriver: mutation   mutation_PDB_numbering   average_ddG   average_total_energy
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
            data.append(d)

    omni[grouping] = data
    aa_order = 'I V L F C M A G T S W Y P H N D E Q K R'.split()
    sdata = {}  # dictionary order is guaranteed with 3.7
    labeldata = {}

    for entry in data:
        if entry['chain'] in (this_chain,):
            if entry['resi'] not in sdata:
                sdata[entry['resi']] = {aa: 0 for aa in aa_order}
                labeldata[entry['resi']] = {}
            sdata[entry['resi']][entry['to_resn']] = entry['ddG']
            labeldata[entry['resi']][entry['to_resn']] = entry['mutation']

    # str_to_normal_form = lambda x: [dict(zip(('from_residue','residue_index','to_residue'), (m[0],m[1:-1],m[-1]))) for m in x.split()]
    # gnomad_svbp = str_to_normal_form('R19K S26* A27S L31V K32R Q35R A37* E38* Y40C V45I V45D M46K E48* L49P E50G Q51* Q53* D55N C58G M61R P63S P64S G65V')
    # gel_svbp = str_to_normal_form('R32R D55N D55A Q28* L31V R32R L49P')

    df = pd.DataFrame.from_records(sdata) \
        .applymap(lambda x: x if x < 20 else 20) \
        .reindex(index=aa_order)

    labels = pd.DataFrame.from_records(labeldata).reindex(index=aa_order)

    omni['labels'] = labels

    muts_to_shapes = lambda muts, color, size_mod=0: [go.layout.Shape(
        type="rect",
        y0=int(v['resi']) - .5 - size_mod,
        x0=aa_order.index(v["to_resn"]) - 0.5 - size_mod,
        y1=int(v['resi']) + .5 + size_mod,
        x1=aa_order.index(v["to_resn"]) + 0.5 + size_mod,
        line=dict(
            color=color,
            width=1
        ),
    ) for v in muts if v["to_resn"] != '*']

    gmuts_to_shapes = lambda muts, color, size_mod=0: [go.layout.Shape(
        type="rect",
        y0=int(v.x) - .5 - size_mod,
        x0=aa_order.index(re.match(r'\w\d+([\w*]+)', v.description).group(1)) - 0.5 - size_mod,
        y1=int(re.match(r'\w(\d+)\w', v.description).group(1)) + .5 + size_mod,
        x1=aa_order.index(re.match(r'\w\d+([\w*]+)', v.description).group(1)) + 0.5 + size_mod,
        line=dict(
            color=color,
            width=1
        ),
    ) for v in muts if v.type == 'missense']

    fig = go.Figure(data=go.Heatmap(x=df.index,
                                    y=df.columns,
                                    z=df.transpose(),
                                    colorscale='temps',
                                    zmid=0,
                                    # reversescale=True,
                                    hovertext=labels))
    fig.update_layout(
        title={'text': title, 'x': 0.5, 'xanchor': 'center'},
        xaxis={'title': 'Residue index'},
        yaxis={'title': 'Residue changed to', 'categoryorder': "trace", 'autorange': 'reversed'},
        height=1600,
        shapes=muts_to_shapes(clinvar['other'], "darkgrey")
               + muts_to_shapes(clinvar['ours'], "black")
               + gmuts_to_shapes(p.gnomAD, "white")
    )
    fig.show()
    fig.write_image('/home/matteo/Desktop/GNB2_part2/' + grouping + '.png', scale=3)