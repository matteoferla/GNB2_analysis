########################################################################################################################
# !/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

if sys.version_info[0] < 3:
    raise EnvironmentError("Hey, caveman, use Python 3.")

__doc__ = \
    """
    Variant plots

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


def model_heatmap():
    df = pd.read_csv('GNB2_analysis.csv')
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

def diff_heatmap():
    df = pd.read_csv('GNB2_analysis.csv')
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
    if mutant in pathogenic:
        return 'pathogenic'
    elif mutant in gnomad:
        return 'gnomad'
    else:
        raise ValueError

def interface_heatmap():
    df = pd.read_csv('GNB2_analysis.csv')
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