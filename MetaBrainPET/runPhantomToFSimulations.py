#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ./runPhantomToFSimulations.py
# The script runs a command as this example: Gate -a [crystal,BGO][acq_time,0.1] mainMacro-PhantoToF.mac > output/PointPhantomToF/BGO/output.txt
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
    materials = ['LYSO-Proteus', 'BGO', 'LYSO-Proteus_BaF2_30', 'BGO_BaF2_25'] #'LYSO-Proteus', 'BGO', 'BGO_BaF2_25', 'BGO_EJ232_25',LYSO-Proteus_EJ232_30
    print(f'PhantomToF source simulations')
    for material in materials:
        print(f'material is {material}')
        time = 120 #in seconds
        time_string = f'{time:.5f}'
        print(f'acquisition time is {time/ 1e0:.5f} s')
        cmd = 'mkdir -p output/PhantomToF/' + material
        output = os.system(cmd)
        print(f'{cmd} ran with exit code {output}')
        cmd = 'Gate -a [crystal,' + material + '][acq_time,' + time_string + '] mainMacro-PhantomToF.mac > ./output/PhantomToF/' + material + '/output.txt &'
        output = os.system(cmd)
        print(f'{cmd} ran with exit code {output}')
    r = True
    return r


# --------------------------------------------------------------------------
if __name__ == '__main__':
    analysis_all_click()
