#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# sudo apt install python3-pip
# pip install uproot
# pip install gatetools
# ./runAnalyticalSensitivityv3.py ./output/Sensitivity > ./output/Sensitivity/AnalyticalSensitivity.txt
# evaluate fitting with just one point (peak sensitivity)
# include Compton+PE scheme
# include off axis sensitivity
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
import math
from array import array
import pandas as pd

DPI = 100

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('output_folder', nargs=-1, required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=True))
@gt.add_options(gt.common_options)

def analysis_all_click(output_folder, **kwargs):
    analyse_all_folders(output_folder)

def analyse_all_folders(output_folder, **kwargs):
    FontSize = 0.05
    ROOT.gStyle.SetLabelSize(FontSize,"xyz")
    ROOT.gStyle.SetTitleSize(FontSize,"xyz")
    ROOT.gStyle.SetTextSize(FontSize)
    ROOT.gStyle.SetLegendTextSize(FontSize)
    ROOT.gStyle.SetLineWidth( 3 )
    ROOT.gStyle.SetLineStyle( 1 )
    ROOT.gStyle.SetMarkerSize( 2 )
    ROOT.gStyle.SetMarkerStyle( 20 )
    ROOT.gROOT.ForceStyle()
    Color = [1,2,3,4,6,9,12,15,28]
    legend_text = [' - AM (511 keV)', ' - AM (^{22}Na)', ' - GATE (^{22}Na)', ' - GATE (511 keV)', ' - AM (only p.e.)', ' - AM (fit #Phi_{scint})']
    BulkArea = 9.0 #mm2
    AreaHZBaF2 = 5*0.3*3 #mm2 entrance area of High Z element combined with BaF2
    AreaBaF2 = 5*0.3*3 #mm2 area of Baf2
    AreaMBaF2 = AreaHZBaF2 + AreaBaF2
    AreaHZEJ232 = 7*0.3*3 #mm2 entrance area of High Z element combined with EJ232
    AreaEJ232 = 8*0.1*3 #mm2
    AreaMEJ232 = AreaHZEJ232 + AreaEJ232
    Scintillators = ['LYSO', 'BGO', 'LYSO-BaF2', 'BGO-BaF2', 'LYSO-EJ232', 'BGO-EJ232']
    n_scint = len(Scintillators)
    Thickness = array( 'f', [ 20., 15., 30., 25., 30., 25. ] )
    EntArea = array( 'f', [ BulkArea, BulkArea, AreaMBaF2, AreaMBaF2, AreaMEJ232, AreaMEJ232 ] )
    AxialPositions = array( 'f', [ -75., -60., -45., -30., -15., 0., 15., 30., 45., 60., 75. ] )
    n_pos = len(AxialPositions)
    LYSO = [("Lu18Y2Si10O50", 7.1, BulkArea, Thickness[0])]
    BGO = [("Bi4Ge3O12", 7.13, BulkArea, Thickness[1])]
    LYSO_BaF2 = [("Lu18Y2Si10O50", 7.1, AreaHZBaF2, Thickness[2]), ("BaF2", 4.88, AreaBaF2, Thickness[2])]
    BGO_BaF2 = [("Bi4Ge3O12", 7.13, AreaHZBaF2, Thickness[3]), ("BaF2", 4.88, AreaBaF2, Thickness[3])]
    LYSO_EJ232 = [("Lu18Y2Si10O50", 7.1, AreaHZEJ232, Thickness[4]), ("C466H513", 1.023, AreaEJ232, Thickness[4])]
    BGO_EJ232 = [("Bi4Ge3O12", 7.13, AreaHZEJ232, Thickness[5]), ("C466H513", 1.023, AreaEJ232, Thickness[5])]
    Material_Scint = [LYSO, BGO , LYSO_BaF2, BGO_BaF2, LYSO_EJ232, BGO_EJ232]
    NRings = 6
    NModPerRing = 30
    NElPerMod = 64
    RingDiam = 290.0 #mm
    AxialExt = 175.2 #mm
    CylArea = 2*math.pi*(RingDiam/2)*AxialExt
    Omega0 = AxialExt/math.sqrt(AxialExt**2+4*(RingDiam/2)**2)
    Omega0Ref = AxialExt/math.sqrt((AxialExt/2.)**2+(RingDiam/2)**2)
    #print('Omega0: '+str(Omega0)+' Omega0Ref: '+str(Omega0Ref))
    index = 0
    start = 0
    stop = 6
    PhiFac, PhiFacErr = array( 'd' ), array( 'd' )
    ref_data_511keVEW = array( 'f')
    ref_data_Na22EW = array( 'f')
    ref_data_511keVnoEW = array( 'f')
    ref_data_Na22noEW = array( 'f')
    pd.options.display.width = 0
    pd.options.display.precision=2
    dfnoEW = pd.DataFrame(columns = ['PeakSnoEW511keVRef', 'PeakSnoEW511keV', 'PeakSnoEW511keVGATE', 'PeakSnoEWNa22', 'PeakSnoEWNa22GATE'], index = ['LYSO', 'BGO', 'LYSO-BaF2', 'BGO-BaF2', 'LYSO-EJ232', 'BGO-EJ232'])
    dfEW = pd.DataFrame(columns = ['Phifac', 'RedChi2', 'PeakSEWPE', 'PeakSEWPhi', 'PeakSEW511keVGATE'], index = ['LYSO', 'BGO', 'LYSO-BaF2', 'BGO-BaF2', 'LYSO-EJ232', 'BGO-EJ232'])
    dfStat = pd.DataFrame(columns = ['MAEnoEWNa22', 'MAPEnoEWNa22', 'STDEnoEWNa22', 'MAEnoEW511keV', 'MPAEnoEW511keV', 'STDEnoEW511keV'], index = ['LYSO', 'BGO', 'LYSO-BaF2', 'BGO-BaF2', 'LYSO-EJ232', 'BGO-EJ232'])
    grAS511keVnoEW=[]
    grASNa22noEW=[]
    grASGATENa22noEW=[]
    grASGATE511keVnoEW=[]
    grAS511keVFitEW=[]
    grAS511keVPheEW=[]
    grASGATENa22EW=[]
    grASGATE511keVEW=[]

    input_filename = 'GateData/Gate511keVnoEW.dat'
    GateS511keVnoEW = read_input(input_filename)
    input_filename = 'GateData/Gate511keVEW.dat'
    GateS511keVEW = read_input(input_filename)
    input_filename = 'GateData/GateNa22noEW.dat'
    GateSNa22noEW = read_input(input_filename)
    input_filename = 'GateData/GateNa22EW.dat'
    GateSNa22EW = read_input(input_filename)

    for j, Scintillator in enumerate(Scintillators[start:stop], start=start):
        print(f' ')
        print('Scintillator: '+Scintillator+' Thick: '+str(Thickness[j]))
        #scint_volume = Thickness[j]*EntArea[j]
        #scint_thicknesses = [10., 12., 15., 18., 20., 22., 25., 28., 30.]
        #dfSPF = pd.DataFrame(columns = ['entrance area', 'packing factor', 'peak sensitivity no EW'])
        #for thick in scint_thicknesses:
        #    scint_ent_area = scint_volume/thick
        #    sens_total_t = pet_performance.sensitivity("PET", Scintillator, round(TotalAttCoef[j],4), thick, RingDiam, AxialExt, NRings, NModPerRing, NElPerMod, scint_ent_area)
        #    if sens_total_t.get_packing_fraction() < 1.: 
        #        dfSPF.loc[thick]=[scint_ent_area, sens_total_t.get_packing_fraction(), sens_total_t.compute_axial("511 keV", 0., 1)]

        print('')
        #print(dfSPF)    
        sens = pet_performance.sensitivity("PET", Material_Scint[j], Thickness[j], RingDiam, AxialExt, NRings, NModPerRing, NElPerMod, EntArea[j])
        ref_data_511keVnoEW = array( 'f', GateS511keVnoEW[j])
        ref_data_Na22noEW = array( 'f', GateSNa22noEW[j])
        ref_data_511keVEW = array( 'f', GateS511keVEW[j])
        ref_data_Na22EW = array( 'f', GateSNa22EW[j])
        ref_data_err = array('f',[0.1 for datai in ref_data_511keVEW])
        S511keVnoEW = array( 'f', sens.compute_axial_array("total", "511 keV", AxialPositions, 1))
        grAS511keVnoEW.append(ROOT.TGraph( n_pos, AxialPositions,  S511keVnoEW))
        SNa22noEW = array( 'f', sens.compute_axial_array("total", "sodium-22", AxialPositions, 1))
        grASNa22noEW.append(ROOT.TGraph( n_pos, AxialPositions, SNa22noEW))
        grASGATENa22noEW.append(ROOT.TGraph( n_pos, AxialPositions, ref_data_Na22noEW))
        grASGATE511keVnoEW.append(ROOT.TGraph( n_pos, AxialPositions, ref_data_511keVnoEW))

        #S511keVPheEW = array( 'f', sens.compute_axial_array("lower_bound","511 keV", AxialPositions, 1))
        S511keVPheEW = array( 'f', sens.compute_axial_array("upper_bound","511 keV", AxialPositions, 1))
        #S511keVPheEW = array( 'f', sens.compute_axial_array("mid_bound","511 keV", AxialPositions, 1))
        #S511keVPheEW = array( 'f', sens.compute_axial_array("photoelectric","511 keV", AxialPositions, 1))

        grAS511keVPheEW.append(ROOT.TGraph( n_pos, AxialPositions, S511keVPheEW))
        grASGATENa22EW.append(ROOT.TGraph( n_pos, AxialPositions, ref_data_Na22EW))
        grASGATE511keVEW.append(ROOT.TGraph( n_pos, AxialPositions, ref_data_511keVEW))
        phi, phierr, reduced_chi2 = sens.fit_axial_curve("total","511 keV", AxialPositions, ref_data_511keVEW, ref_data_err)
        S511keVFitEW = array( 'f', sens.compute_axial_array("total","511 keV", AxialPositions, phi))
        grAS511keVFitEW.append(ROOT.TGraph( n_pos, AxialPositions, S511keVFitEW))

        dfSProfiles = pd.DataFrame(columns = ['SnoEW511keV', 'SnoEW511keVGATE', 'SnoEWNa22', 'SnoEWNa22GATE', 'SEW511keVPhe', 'SEW511keVFit', 'SEW511keVGATE', 'SEWNa22GATE'])
        for i, AxialPos in enumerate(AxialPositions, start=0): dfSProfiles.loc[AxialPos] = [S511keVnoEW[i], ref_data_511keVnoEW[i], SNa22noEW[i], ref_data_Na22noEW[i], S511keVPheEW[i], S511keVFitEW[i], ref_data_511keVEW[i], ref_data_Na22EW[i]]
        print(dfSProfiles)
        dfnoEW.loc[Scintillator] = [sens.get_peak_sensitivity_ref(), sens.compute_axial("total","511 keV", 0., 1), max(ref_data_511keVnoEW), sens.compute_axial("total","sodium-22", 0., 1), max(ref_data_Na22noEW)]
        dfEW.loc[Scintillator] = [phi, reduced_chi2, sens.compute_axial("upper_bound","511 keV", 0., 1), sens.compute_axial("total","511 keV", 0., phi), max(ref_data_511keVEW)]
        
        STDnoEW511keV = np.std(np.subtract(ref_data_511keVnoEW,S511keVnoEW))
        STDnoEWNa22 = np.std(np.subtract(ref_data_Na22noEW,SNa22noEW))
        dfStat.loc[Scintillator]=[mae(ref_data_Na22noEW,SNa22noEW),mape(ref_data_Na22noEW,SNa22noEW),STDnoEWNa22,mae(ref_data_511keVnoEW,S511keVnoEW),mape(ref_data_511keVnoEW,S511keVnoEW),STDnoEW511keV]

    print('')
    print(dfnoEW)
    print(dfEW)
    print(dfStat)
    c1 = ROOT.TCanvas( 'c1', '', 1700, 1200 )
    c1.SetGridx(1)
    c1.SetGridy(0)
    c1.SetRightMargin(0.02)
    c1.SetTopMargin(0.08)
    c1.SetBottomMargin(0.15)
    c1.SetLeftMargin(0.12)

    for k in range(start, stop, 2):
        MgrASnoEW = ROOT.TMultiGraph()
        MgrASnoEW.SetTitle(';axial position (mm);sensitivity (%)')
        MgrASnoEW.SetMaximum(80.)
        SetLineAndMarkerColorAndLegend(grAS511keVnoEW[k], Color[0], Scintillators[k]+legend_text[0])
        MgrASnoEW.Add(grAS511keVnoEW[k],'L')
        SetLineAndMarkerColorAndLegend(grAS511keVnoEW[k+1], Color[2], Scintillators[k+1]+legend_text[0])
        MgrASnoEW.Add(grAS511keVnoEW[k+1],'L')
        SetLineAndMarkerColorAndLegend(grASGATE511keVnoEW[k], Color[0],Scintillators[k]+legend_text[3])
        MgrASnoEW.Add(grASGATE511keVnoEW[k],'P')
        SetLineAndMarkerColorAndLegend(grASGATE511keVnoEW[k+1], Color[2],Scintillators[k+1]+legend_text[3])
        MgrASnoEW.Add(grASGATE511keVnoEW[k+1],'P')
        SetLineAndMarkerColorAndLegend(grASNa22noEW[k], Color[1], Scintillators[k]+legend_text[1])
        MgrASnoEW.Add(grASNa22noEW[k],'L')
        SetLineAndMarkerColorAndLegend(grASNa22noEW[k+1], Color[3], Scintillators[k+1]+legend_text[1])
        MgrASnoEW.Add(grASNa22noEW[k+1],'L')
        SetLineAndMarkerColorAndLegend(grASGATENa22noEW[k], Color[1],Scintillators[k]+legend_text[2])
        MgrASnoEW.Add(grASGATENa22noEW[k],'P')
        SetLineAndMarkerColorAndLegend(grASGATENa22noEW[k+1], Color[3],Scintillators[k+1]+legend_text[2])
        MgrASnoEW.Add(grASGATENa22noEW[k+1],'P')
        MgrASnoEW.Draw('A')
        legnoEW = ROOT.TLegend(0.44, 0.55, 0.99, 0.99)
        legnoEW.AddEntry(grASGATE511keVnoEW[k],'','P')
        legnoEW.AddEntry(grASGATENa22noEW[k],'','P')
        legnoEW.AddEntry(grAS511keVnoEW[k],'','L')
        legnoEW.AddEntry(grASNa22noEW[k],'','L')
        legnoEW.AddEntry(grASGATE511keVnoEW[k+1],'','P')
        legnoEW.AddEntry(grASGATENa22noEW[k+1],'','P')
        legnoEW.AddEntry(grAS511keVnoEW[k+1],'','L')
        legnoEW.AddEntry(grASNa22noEW[k+1],'','L')
        legnoEW.Draw('same')
        c1.Update()
        c1.SaveAs(output_folder[0]+'/sens-'+Scintillators[k]+'-'+Scintillators[k+1]+'noEW.pdf')

        MgrASEW = ROOT.TMultiGraph()
        MgrASEW.SetTitle(';axial position (mm);sensitivity (%)')
        MgrASEW.SetMaximum(20.)
        SetLineAndMarkerColorAndLegend(grASGATE511keVEW[k], Color[0], Scintillators[k]+legend_text[3])
        MgrASEW.Add(grASGATE511keVEW[k],'P')
        SetLineAndMarkerColorAndLegend(grAS511keVPheEW[k], Color[0], Scintillators[k]+legend_text[4])
        MgrASEW.Add(grAS511keVPheEW[k],'L')
        SetLineAndMarkerColorAndLegend(grAS511keVFitEW[k], Color[0], Scintillators[k]+legend_text[5])
        grAS511keVFitEW[k].SetLineStyle(9)
        MgrASEW.Add(grAS511keVFitEW[k],'L')
        SetLineAndMarkerColorAndLegend(grASGATENa22EW[k], Color[1], Scintillators[k]+legend_text[2])
        MgrASEW.Add(grASGATENa22EW[k],'P')
        SetLineAndMarkerColorAndLegend(grASGATENa22EW[k+1], Color[3], Scintillators[k+1]+legend_text[2])
        MgrASEW.Add(grASGATENa22EW[k+1],'P')
        SetLineAndMarkerColorAndLegend(grASGATE511keVEW[k+1], Color[2], Scintillators[k+1]+legend_text[3])
        MgrASEW.Add(grASGATE511keVEW[k+1],'P')
        SetLineAndMarkerColorAndLegend(grAS511keVFitEW[k+1], Color[2], Scintillators[k+1]+legend_text[5])
        grAS511keVFitEW[k+1].SetLineStyle(9)
        MgrASEW.Add(grAS511keVFitEW[k+1],'L')
        SetLineAndMarkerColorAndLegend(grAS511keVPheEW[k+1], Color[2], Scintillators[k+1]+legend_text[4])
        MgrASEW.Add(grAS511keVPheEW[k+1],'L')
        MgrASEW.Draw('A') 
        legEW = ROOT.TLegend(0.45, 0.55, 0.99, 0.99)
        legEW.AddEntry(grASGATE511keVEW[k],'','P')
        legEW.AddEntry(grASGATENa22EW[k],'','P')
        legEW.AddEntry(grAS511keVPheEW[k],'','L')
        legEW.AddEntry(grAS511keVFitEW[k],'','L')
        legEW.AddEntry(grASGATE511keVEW[k+1],'','P')
        legEW.AddEntry(grASGATENa22EW[k+1],'','P')
        legEW.AddEntry(grAS511keVPheEW[k+1],'','L')
        legEW.AddEntry(grAS511keVFitEW[k+1],'','L')
        legEW.Draw('same')
        c1.Update()
        c1.SaveAs(output_folder[0]+'/sens-'+Scintillators[k]+'-'+Scintillators[k+1]+'EW.pdf')


    pde = pet_performance.photopeak_efficiency(LYSO, 20.)
    # A grid of scattering angles in rad.
    theta = np.arange(0, np.pi, 0.01)
    n = len(theta)
    KN_dsigma_dTheta = pde.get_KN_dsigma_dTheta(theta, 511e3)
    grKN_dsigma_dTheta = ROOT.TGraph(n, np.degrees(theta), KN_dsigma_dTheta)
    grKN_dsigma_dTheta.SetTitle(';theta (deg); dsigma_dTheta')
    grKN_dsigma_dTheta.Draw('ACP')
    c1.Update()
    c1.SaveAs(output_folder[0]+'/KN_dsigma_dTheta.pdf')

    KN_dsigma_dOmega = pde.get_KN_dsigma_dOmega(theta, 511e3)
    grKN_dsigma_dOmega = ROOT.TGraph(n, np.degrees(theta), KN_dsigma_dOmega)
    grKN_dsigma_dOmega.SetTitle(';theta (deg); dsigma_dOmega')
    grKN_dsigma_dOmega.Draw('ACP')
    c1.Update()
    c1.SaveAs(output_folder[0]+'/KN_dsigma_dOmega.pdf')


    r = True
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

def read_input(filename):
    data = []
    with open(filename, 'r') as file:
        for line in file:
            # Remove comments at the end of each line
            line = line.split('#')[0].strip()
            if line:  # Ignore empty lines
                data.append(list(map(float, line.split(','))))
    return np.array(data)

# --------------------------------------------------------------------------
if __name__ == '__main__':
    analysis_all_click()

