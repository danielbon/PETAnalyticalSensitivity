#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ./runPointSourceSimulations.py
# The script runs a command as this example: Gate -a [crystal,BGO][axial_pos,0][acq_time,0.1] mainMacro-PointSource.mac > output/PointSource/BGO/0/output.txt
from __future__ import print_function

import click
import os
import numpy as np
#import gatetools as gt

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('output_folders', nargs=-1, required=False,
                type=click.Path(exists=True, file_okay=True, dir_okay=True))
#@gt.add_options(gt.common_options)
def analysis_all_click(output_folders, **kwargs):
    analyse_all_folders(output_folders)

def analyse_all_folders(output_folders, **kwargs):
    axial_pos = ['0'] #, '5'] # in mm
    materials = ['LYSO-Proteus', 'BGO', 'BGO_BaF2_25', 'BGO_EJ232_25']
    n_axial_pos = len(axial_pos)
    print(f'Point source simulations')
    for material in materials:
        print(f'material is {material}')
        for position in axial_pos:
            print(f'position is {position}')
            time = 10 #in seconds
            time_string = f'{time:.5f}'
            print(f'acquisition time is {time/ 1e0:.5f} s')
            cmd = 'mkdir -p output/PointSource/' + material + '/' + position
            output = os.system(cmd)
            print(f'{cmd} ran with exit code {output}')
            cmd = 'Gate -a [crystal,' + material + '][axial_pos,' + position + '][acq_time,' + time_string + '] mainMacro-PointSource.mac > ./output/PointSource/' + material + '/' + position + '/output.txt &'
            output = os.system(cmd)
            print(f'{cmd} ran with exit code {output}')
            
    r = True
    return r


# --------------------------------------------------------------------------
if __name__ == '__main__':
    analysis_all_click()
