########################################################################################################################
# !/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

if sys.version_info[0] < 3:
    raise EnvironmentError("Hey, caveman, use Python 3.")

__doc__ = \
    """
    Variant plots... these are meant to be run in a Jupyter notebook and are here just because.

    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2019 A.D."
__license__ = "MIT"
__version__ = "1"
__citation__ = "TBA"

########################################################################################################################

import pandas as pd
import plotly.graph_objects as go


def model_heatmap(datafile:str):
    df = pd.read_csv(datafile)
    for n in ('pathogenic', 'gnomad'):
        dfixed = df.loc[df.groupname == n]
        ddG = dfixed.pivot_table(index=['model'], columns=['mutation'], values='ddG').round(1)

        fig = go.Figure(data=go.Heatmap(
            z=ddG,
            x=ddG.columns,
            y=ddG.index,
            hoverongaps=False,
            colorscale='temps'))
        fig.data[0].update(zmin=-20, zmax=20)
        fig.update_layout(title=f'ΔΔG relative to model used for {n} variants')
        fig.write_image('heatmap_' + n + '.png', scale=3)
        fig.show()

def diff_heatmap(datafile:str):
    df = pd.read_csv(datafile)
    for n in ('pathogenic', 'gnomad'):
        dfixed = df.loc[df.groupname == n]
        ddG = dfixed.pivot_table(index=['mutation'], columns=['model'], values='ddG').round(1)
        ddG = ddG.apply(lambda row: row.apply(lambda n: n - row.unpaired), 1)
        # ddG = ddG.assign(ddG=ddG.ddG.apply(lambda n: -20 if n < -20 else n).apply(lambda n: 20 if n > 20 else n))

        # new_order = list(sorted(ddG.index.values, key=lambda n: int(n[1:-1])))
        ddG = ddG[['wAlpha', 'wGRK2', 'wPREX1', 'phospho', 'wKCTD12']]
        # .sort_values(new_order,axis=1)
        print('ddG')
        display(ddG)

        fig1 = go.Figure(data=go.Heatmap(
            z=ddG.transpose(),
            x=ddG.index,
            y=ddG.columns,
            hoverongaps=False,
            colorscale='temps'))
        fig1.data[0].update(zmin=-20, zmax=20)
        fig1.update_layout(
            title='Difference in free energy (ΔΔG, REU, approx. kcal/mol) for variants relative to unbound form')
        fig1.write_image('heatmap_crudeDiff_' + n + '.png', scale=3)
        fig1.show()

def get_group(mutant):
    gnomad = ['A21S', 'R49K', 'A73T', 'R129H', 'V133I', 'D195N', 'R197H', 'T198M', 'I208V', 'V213M', 'M263V', 'R283W',
              'G302S', 'D303N', 'A309T', 'D322H', 'D323N', 'S2R', 'E12D', 'R22Q', 'T31I', 'G36R', 'G36E', 'D38E',
              'I43M', 'R49K', 'D66E', 'V71L', 'L79V', 'S84T', 'N88T', 'A104T', 'Y105F', 'C114S', 'I120V', 'I123V',
              'I123T', 'R129H', 'R137G', 'T143A', 'I157V', 'T159S', 'T164A', 'T173I', 'G174S', 'V178L', 'S191A',
              'A193S', 'D195N', 'R197C', 'R197H', 'T198K', 'S207C', 'I208V', 'K209R', 'D212H', 'V213M', 'M217V',
              'R219Q', 'I223V', 'F241V', 'G244V', 'A248V', 'T249M', 'F253L', 'D258Y', 'L262V', 'M263V', 'M263T',
              'H266N', 'H266R', 'N268D', 'N268I', 'G272S', 'S275A', 'R280C', 'R280H', 'A287T', 'I296T', 'A299T',
              'M300V', 'G302S', 'D303N', 'R304H', 'A305T', 'A309G', 'D312V', 'V315M', 'L318I', 'G319X', 'D322N',
              'D322G', 'D323N', 'M325V', 'V327M', 'F335L', 'I338V']
    pathogenic = ('R52L', 'A73T', 'G77R', 'G77E', 'K89T', 'K89E', 'E180K', 'S147L', 'I171T')
    if mutant in pathogenic:
        return 'pathogenic'
    elif mutant in gnomad:
        return 'gnomad'
    else:
        raise ValueError

def interface_heatmap(datafile:str):
    df = pd.read_csv(datafile)
    ref = {row.mutation: row.interface_Δscore for i, row in df.iterrows() if row.model == 'unpaired'}
    # display(df.loc[df.mutation == 'S147L'])
    dfixed = df.assign(score=df.apply(lambda row: row.interface_Δscore - ref[row.mutation], 1)) \
        .loc[df.groupname == 'pathogenic'] \
        .loc[df.model != 'unpaired']
    # dfixed = dfixed.assign(score=dfixed.score.apply(lambda n: -20 if n < -20 else n).apply(lambda n: 20 if n > 20 else n))
    interface = dfixed.pivot_table(index=['model'], columns=['mutation'], values='score').round(1)
    print('interface')
    display(interface)

    fig = go.Figure(data=go.Heatmap(
        z=interface,
        x=interface.columns,
        y=interface.index,
        hoverongaps=False,
        colorscale='temps'))
    fig.data[0].update(zmin=-20, zmax=20)
    fig.update_layout(
        title='Difference in interface energy (ΔΔG, REU, approx. kcal/mol) for variants relative to unbound form')
    fig.write_image('heatmap_interface.png', scale=3)
    fig.show()