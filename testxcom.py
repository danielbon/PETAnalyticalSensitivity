#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# sudo apt install python3-pip
# pip install uproot
# pip install gatetools
# ./runAnalyticalSensitivityv2.py ./output/Sensitivity > ./output/Sensitivity/AnalyticalSensitivity.txt
from __future__ import print_function

import pet_performance
import uproot
import click
import gatetools as gt
import os
import numpy as np
from box import Box
from statistics import mean
import ROOT
from ROOT import TMath
import math
from math import sin
from array import array
import pandas as pd
from nistcalculators.xcom import xcom


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('output_folder', nargs=-1, required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=True))
@gt.add_options(gt.common_options)

def analysis_all_click(output_folder, **kwargs):
    analyse_all_folders(output_folder)

def analyse_all_folders(output_folder, **kwargs):

    BulkArea = 9.0 #mm2
    AreaHZBaF2 = 5*0.3*3 #mm2 entrance area of High Z element combined with BaF2
    AreaBaF2 = 5*0.3*3 #mm2 area of Baf2
    AreaMBaF2 = AreaHZBaF2 + AreaBaF2
    AreaHZEJ232 = 7*0.3*3 #mm2 entrance area of High Z element combined with EJ232
    AreaEJ232 = 8*0.1*3 #mm2
    AreaMEJ232 = AreaHZEJ232 + AreaEJ232
    Scintillators = ['LYSO', 'BGO', 'LYSO-BaF2', 'BGO-BaF2', 'LYSO-EJ232', 'BGO-EJ232']
    n_scint = len(Scintillators)
    TotalAttCoef = array( 'f', [ 0.819, 0.963, 0.638, 0.710, 0.619, 0.723 ] )
    PhEAbsCoef = array( 'f', [ 0.254, 0.396, 0.169, 0.241, 0.184, 0.287 ] )
    Thickness = array( 'f', [ 20., 15., 30., 25., 30., 25. ] )
    EntArea = array( 'f', [ BulkArea, BulkArea, AreaMBaF2, AreaMBaF2, AreaMEJ232, AreaMEJ232 ] )

    LYSO = ["Lu18Y2Si10O50", 7.1, BulkArea, Thickness[0]]
    BGO = ["Bi4Ge3O12", 7.1, BulkArea, Thickness[1]]
    LYSO_BaF2 = [("Lu18Y2Si10O50", 7.1, AreaHZBaF2, Thickness[2]), ("BaF2", 4.88, AreaBaF2, Thickness[2])]
    BGO_BaF2 = [("Bi4Ge3O12", 7.1, AreaHZBaF2, Thickness[3]), ("BaF2", 4.88, AreaBaF2, Thickness[3])]
    LYSO_EJ232 = [("Lu18Y2Si10O50", 7.1, AreaHZEJ232, Thickness[4]), ("C466H513", 1.023, AreaEJ232, Thickness[4])]
    BGO_EJ232 = [("Bi4Ge3O12", 7.1, AreaHZEJ232, Thickness[5]), ("C466H513", 1.023, AreaEJ232, Thickness[5])]

    #sens_total = pet_performance.sensitivity("PET", Scintillator, Thickness[j], RingDiam, AxialExt, NRings, NModPerRing, NElPerMod, EntArea[j])
    m_total = 0.
    area_total = 0.
    for i in range(len(LYSO_BaF2)): 
        m_total += LYSO_BaF2[i][1]*LYSO_BaF2[i][2]*LYSO_BaF2[i][3]
        area_total += LYSO_BaF2[i][2]
    print(m_total, area_total)
    formulas_scint=[]
    weights_scint=[]
    for i in range(len(LYSO_BaF2)): 
        weights_scint.append((LYSO_BaF2[i][1]*LYSO_BaF2[i][2]*LYSO_BaF2[i][3])/m_total)
        formulas_scint.append(LYSO_BaF2[i][0])
    eff_density = m_total/(area_total*Thickness[2])
    matLYSO_BaF2 = xcom.MaterialFactory.mix_formulas(formulas=formulas_scint, weights=weights_scint)
    data_LYSO_BaF2 = xcom.calculate_attenuation(matLYSO_BaF2, [511e3, 1274e3])
    print(eff_density)
    print(matLYSO_BaF2.elements_by_Z, matLYSO_BaF2.weights)
    print(data_LYSO_BaF2)
    print(data_LYSO_BaF2[0]['total']*eff_density)
    print(data_LYSO_BaF2[0]['photoelectric']*eff_density)

    densLYSO = 7.1
    matLYSO = xcom.MaterialFactory.from_formula("Lu18Y2Si10O50")
    data_LYSO = xcom.calculate_attenuation(matLYSO, [511e3, 1274e3])
    print(matLYSO.elements_by_Z, matLYSO.weights)
    print(data_LYSO)
    for name in data_LYSO.dtype.names:
        if name != "energy": data_LYSO[name] *= (densLYSO)
    print(data_LYSO['incoherent'])
    densBGO = 7.1
    BGO = xcom.MaterialFactory.from_formula("Bi4Ge3O12")
    data_BGO = xcom.calculate_attenuation(BGO, [511e3, 1274e3])
    print(BGO.elements_by_Z, BGO.weights)
    print(data_BGO)
    for name in data_BGO.dtype.names:
        if name != "energy": data_BGO[name] *= (densBGO)
    print(data_BGO)
    print(data_BGO[1]['incoherent'])
    densBaF2 = 4.88
    BaF2 = xcom.MaterialFactory.from_formula("BaF2")
    data_BaF2 = xcom.calculate_attenuation(BaF2, [511e3, 1274e3])
    print(BaF2.elements_by_Z, BaF2.weights)
    print(data_BaF2)
    for name in data_BaF2.dtype.names:
        if name != "energy": data_BaF2[name] *= (densBaF2)
    print(data_BaF2)
    print(data_BaF2[1]['incoherent'])

    #xcom.calculate_attenuation(xcom.MaterialFactory.from_formula("H2O"), [1e3])
    #array([(1000., 1.3722503, 0.01319334, 4075.62851585, 0., 0., 4075.64170918, 4077.01395948)],
    #dtype=[('energy', '<f8'), ('coherent', '<f8'), ('incoherent', '<f8'), ('photoelectric', '<f8'), ('pair_atom', '<f8'), ('pair_electron', '<f8'), ('total_without_coherent', '<f8'), ('total', '<f8')])


    r = 1
    return r

#Mean absolute percentage error
def mape(y_true, predictions):
    y_true, predictions = np.array(y_true), np.array(predictions)
    return 100.*np.mean(np.abs((y_true - predictions)/y_true))

#Mean absolute error
def mae(y_true, predictions):
    y_true, predictions = np.array(y_true), np.array(predictions)
    return np.mean(np.abs(y_true - predictions))

def SetLineAndMarkerColorAndLegend(graph, color, legend):
    graph.SetLineColor( color )
    graph.SetMarkerColor( color )
    graph.SetTitle(legend)
    return 0

# --------------------------------------------------------------------------
if __name__ == '__main__':
    analysis_all_click()

