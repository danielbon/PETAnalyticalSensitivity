#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ./runNECRSimulations.py
# The script runs a command as this example: Gate -a [crystal,$Material][activity,$activity][acq_time,0.005] mainMacro-NECR.mac&
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
    #activity_values = [ '3000', '2000', '1000', '500', '300', '100', '50', '10'] # in MBq
    activity_values = [ '2500', '1500', '750', '400'] # in MBq
    materials = ['LYSO-Proteus', 'BGO',  'LYSO-Proteus_BaF2_31.1', 'LYSO-Proteus_EJ232_32.3', 'BGO_BaF2_23.7', 'BGO_EJ232_22.7']#
    n_activities = len(activity_values)
    half_life = 6586  # in s
    print(f'NECR simulations')
    print(f'Half life         = {half_life / 1e0:.2f} s ###')
    decay_const = np.log(2)/half_life
    det_evts = 2e6 # number of prompt events to be registered
    print(f'detected events   = {det_evts / 1e0:f} events')
    sensitivity = 0.0149
    print(f'sensitivity       = {sensitivity / 1e0:f} ')
    for material in materials:
        print(f'material is {material}')
        for activity in activity_values:
            time = -np.log(1-det_evts*decay_const/(float(activity)*1e6*sensitivity))/decay_const
            time_string = f'{time:.5f}'
            print(f'activity is {activity}')
            print(f'acquisition time is {time/ 1e0:.5f} s')
            cmd = 'mkdir -p output/NECR/' + material + '/' + activity
            output = os.system(cmd)
            print(f'{cmd} ran with exit code {output}')
            cmd = 'Gate -a [crystal,' + material + '][activity,' + activity + '][acq_time,' + time_string + '] mainMacro-NECR.mac > ./output/NECR/' + material + '/' + activity + '/output.txt &'
            output = os.system(cmd)
            print(f'{cmd} ran with exit code {output}')
            
    r = True
    return r


# --------------------------------------------------------------------------
if __name__ == '__main__':
    analysis_all_click()
