#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# sudo apt install python3-pip
# pip install uproot
# pip install gatetools
# ./runAnalyticalSensitivityv1.1.py ./output/Sensitivity > ./output/Sensitivity/AnalyticalSensitivity.txt
from __future__ import print_function

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
    c1 = ROOT.TCanvas( 'c1', '', 1700, 1200 )
    c1.SetGridx(1)
    c1.SetGridy(0)
    c1.SetRightMargin(0.02)
    c1.SetTopMargin(0.08)
    c1.SetBottomMargin(0.15)
    c1.SetLeftMargin(0.12)
    Color = [1,2,3,4,6,9,12,15,28]
    legend_text = [' - AM (511 keV)', ' - AM (^{22}Na)', ' - GATE (^{22}Na)', ' - GATE (511 keV)', ' - AM (only p.e.)', ' - AM (fit #Phi_{scint})']
    GateS22NanoEW = [[ 12.33, 17.33, 22.23, 26.74, 30.81, 33.72, 30.80, 26.71, 22.22, 17.31, 12.32], #LYSO
                     [ 9.90, 14.03, 18.05, 21.81, 25.21, 27.83, 25.23, 21.80, 18.06, 14.01, 9.90 ], #BGO
                     [ 15.36, 21.63, 27.70, 33.24, 38.24, 41.50, 38.22, 33.23, 27.68, 21.61, 15.37 ], #LYSO_BaF2 30
                     [ 13.31, 18.82, 24.19, 29.10, 33.60, 36.63, 33.51, 29.07, 24.18, 18.78, 13.31 ], #BGO_BaF2 25
                     [ 14.09, 19.87, 25.64, 30.78, 35.59, 38.77, 35.54, 30.79, 25.59, 19.85, 14.10 ], #LYSO_EJ232 30
                     [ 12.54, 17.87, 23.08, 27.95, 32.34, 35.44, 32.38, 27.95, 23.10, 17.89, 12.55 ]] #BGO_EJ232 25
    GateS22NaEW = [[ 2.05, 3.46, 4.85, 6.20, 7.45, 8.47, 7.45,   6.19,   4.85,  3.46,  2.06],#LYSO
                   [ 1.97, 3.41, 4.85, 6.27, 7.60, 8.76, 7.61, 6.27, 4.85, 3.40, 1.97 ],#BGO
                   [ 2.20, 3.62, 5.00, 6.33, 7.56, 8.47, 7.53, 6.31, 4.99, 3.63, 2.21 ],#LYSO_BaF2 30
                   [ 2.27, 3.84, 5.37, 6.85, 8.23, 9.29, 8.20, 6.84, 5.36, 3.83, 2.27 ],#BGO_BaF2 25
                   [ 2.08, 3.42, 4.72, 6.01, 7.20, 8.13, 7.20, 6.01, 4.72, 3.43, 2.08 ],#LYSO_EJ232 30
                   [ 2.32, 3.93, 5.51, 7.10, 8.59, 9.78, 8.59, 7.10, 5.52, 3.96, 2.31 ]]#BGO_EJ232 25
    GateS511keVnoEW = [[ 4.04, 8.16, 12.24, 16.08, 19.82, 22.76, 19.83, 16.10, 12.23, 8.14, 4.03 ],#LYSO
                       [ 3.54, 7.01, 10.44, 13.69, 16.88, 19.58, 16.88, 13.72, 10.44, 6.99, 3.53 ],#BGO
                       [ 4.74, 9.80, 14.82, 19.57, 24.06, 27.31, 24.02, 19.52, 14.84, 9.78, 4.75 ],#LYSO_BaF2 30
                       [ 4.31, 8.83, 13.31, 17.52, 21.62, 24.62, 21.53, 17.49, 13.29, 8.80, 4.29 ],#BGO_BaF2 25
                       [ 4.42, 9.06, 13.80, 18.29, 22.61, 25.78, 22.63, 18.29, 13.81, 9.04, 4.43 ],#LYSO_EJ232 30
                       [ 4.20, 8.50, 12.87, 17.09, 21.13, 24.24, 21.15, 17.09, 12.90, 8.48, 4.19 ] ]#BGO_EJ232 25
    GateS511keVEW = [[ 1.15, 2.49, 3.82, 5.12, 6.37, 7.48, 6.38, 5.12, 3.82, 2.49, 1.15 ],#LYSO
                     [ 1.23, 2.64, 4.05, 5.44, 6.81, 8.07, 6.81, 5.45, 4.06, 2.63, 1.23 ],#BGO
                     [ 1.10, 2.43, 3.72, 4.99, 6.18, 7.13, 6.15, 4.97, 3.71, 2.42, 1.11 ],#LYSO_BaF2 30
                     [ 1.24, 2.74, 4.18, 5.62, 6.99, 8.10, 6.95, 5.60, 4.17, 2.72, 1.24 ],#BGO_BaF2 25
                     [ 1.10, 2.29, 3.51, 4.73, 5.92, 6.90, 5.92, 4.73, 3.50, 2.29, 1.10 ],#LYSO_EJ232 30
                     [ 1.37, 2.85, 4.34, 5.89, 7.37, 8.67, 7.38, 5.89, 4.35, 2.85, 1.38 ]]#BGO_EJ232 25
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
    AxialPositions = array( 'f', [ -75., -60., -45., -30., -15., 0., 15., 30., 45., 60., 75. ] )
    n_pos = len(AxialPositions)
    NRings = 6
    NModPerRing = 30
    NElPerMod = 64
    RingDiam = 290.0 #mm
    AxialExt = 175.2 #mm
    CylArea = 2*math.pi*(RingDiam/2)*AxialExt
    Omega0 = AxialExt/math.sqrt(AxialExt**2+4*(RingDiam/2)**2)
    Omega0Ref = AxialExt/math.sqrt((AxialExt/2.)**2+(RingDiam/2)**2)
    #print('Omega0: '+str(Omega0)+' Omega0Ref: '+str(Omega0Ref))
    BF22Na = 0.9
    index = 0
    start = 0
    stop = 6
    PSensnoEW511keVRef, PSensnoEW511keV, PSensnoEW511keVGATE, PSensnoEW22Na, PSensnoEW22NaGATE = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )
    PSensEWPhi, PSensEWPE, PSensEW511keVGATE, PSensEW22NaGATE = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )
    RedChi2 = array( 'd' )
    PhiFac, PhiFacErr = array( 'd' ), array( 'd' )
    ref_data_511keVEW = array( 'f')
    ref_data_22NaEW = array( 'f')
    ref_data_511keVnoEW = array( 'f')
    ref_data_22NanoEW = array( 'f')
    MAEnoEW511keV, MAEnoEW22Na = array( 'd' ), array( 'd' )
    MAPEnoEW511keV, MAPEnoEW22Na = array( 'd' ), array( 'd' )
    STDEnoEW511keV, STDEnoEW22Na = array( 'd' ), array( 'd' )

    for j, Scintillator in enumerate(Scintillators[start:stop], start=start):
        PF = NRings*NModPerRing*NElPerMod*EntArea[j]/CylArea
        PFRef = NRings*NModPerRing*NElPerMod*EntArea[j]*Thickness[j]/(AxialExt*math.pi*((RingDiam/2+Thickness[j])**2-(RingDiam/2)**2))
        Eipe = 1 - math.exp(-PhEAbsCoef[j]*Thickness[j]*0.1)
        Eitotal = 1 - math.exp(-TotalAttCoef[j]*Thickness[j]*0.1)
        print(f' ')
        print('Scintillator: '+Scintillator+' Thick: '+str(Thickness[j])+' PF: '+str(round(PF,3))+' Eipe: '+str(round(Eipe,3))+' Eitotal: '+str(round(Eitotal,3)))
        pos, SensCorrnoEW = array( 'd' ), array( 'd' )
        posEW, SensCorrEW = array( 'd' ), array( 'd' )
        pos22Na, Sens22NanoEW = array( 'd' ), array( 'd' )
        for i, AxialPos in enumerate(AxialPositions, start=0):
            OmegaAxUncor = 0.5*(((AxialExt/2.)-AxialPos)/math.sqrt(((AxialExt/2.)-AxialPos)**2+(RingDiam/2.)**2)+((AxialExt/2.)+AxialPos)/math.sqrt(((AxialExt/2.)+AxialPos)**2+(RingDiam/2.)**2))          
            OmegaAxCor = ((AxialExt/2.)-abs(AxialPos))/math.sqrt(((AxialExt/2.)-abs(AxialPos))**2+(RingDiam/2.)**2)
            #print('OmegaAxCor: ' + str(round(OmegaAxCor,4)) + '   AxialPos: ' + str(AxialPos))
            SCorrTotnoEW = 100*OmegaAxCor*PF*Eitotal**2
            SensCorrnoEW.append( SCorrTotnoEW )
            S22NanoEW = BF22Na*100*(OmegaAxCor*PF*Eitotal**2+2*(OmegaAxUncor*PF)**2*Eitotal**2)
            Sens22NanoEW.append( S22NanoEW )
            SCorrPEEW = 100*OmegaAxCor*PF*Eipe**2
            SensCorrEW.append( SCorrPEEW )
            posEW.append( AxialPos )
            pos.append( AxialPos )
            if (AxialPos==0.):
                PSensnoEW511keVRef.append(100*Omega0Ref*PFRef**2*Eitotal**2)
                PSensnoEW511keV.append(SCorrTotnoEW)
                PSensnoEW511keVGATE.append( max(array( 'f', GateS511keVnoEW[j])))
                PSensnoEW22Na.append(S22NanoEW)
                PSensnoEW22NaGATE.append(max(array( 'f', GateS22NanoEW[j])))
                PSensEWPE.append(SCorrPEEW)
                PSensEW511keVGATE.append(max(array( 'f', GateS511keVEW[j])))
                PSensEW22NaGATE.append(max(array( 'f', GateS22NaEW[j])))
        print(f'Position(mm)  Sensitivity for correlated events (511 keV) without energy window') 
        for i in range( n_pos ):
            print(' %.1f   %.2f ' % (pos[i],SensCorrnoEW[i]))
        print(f'Position(mm)  Sensitivity for 22Na without energy window') 
        for i in range( n_pos ):
            print(' %.1f   %.2f ' % (pos[i],Sens22NanoEW[i]))
        print(f'Position(mm)  Sensitivity for correlated events (511 keV) with energy window') 
        for i in range( n_pos ):
            print(' %.1f   %.2f ' % (posEW[i],SensCorrEW[i]))

        MgrAxialSensEW = ROOT.TMultiGraph()
        Title = Scintillator + '(EW) ;axial position (mm);sensitivity (%)'
        MgrAxialSensEW.SetTitle(Title)
        MgrAxialSensEW.SetMaximum(15.)
        grASensEW = ROOT.TGraph( n_pos, posEW, SensCorrEW )
        grASensEW.SetLineColor( Color[4] )
        grASensEW.SetMarkerColor( Color[4] )
        grASensEW.SetTitle(legend_text[4])
        GateData = array( 'f', GateS22NaEW[j])
        grASensGATE22NaEW = ROOT.TGraph( n_pos, AxialPositions, GateData )
        grASensGATE22NaEW.SetLineColor( Color[1] )
        grASensGATE22NaEW.SetMarkerColor( Color[1] )
        grASensGATE22NaEW.SetTitle(legend_text[2])
        GateData = array( 'f', GateS511keVEW[j])
        GateDataErr = array('f',[0.1 for GateDatai in GateData])
        grASensGATE511keVEW = ROOT.TGraphErrors( n_pos, AxialPositions, GateData, 0, GateDataErr )
        grASensGATE511keVEW.SetLineColor( Color[0] )
        grASensGATE511keVEW.SetMarkerColor( Color[0] )
        grASensGATE511keVEW.SetTitle(legend_text[3])
        fSensCorrEW = ROOT.TF1("fSensCorrEW", "100*(([0]/2.)-TMath::Abs(x))/TMath::Sqrt((([0]/2.)-TMath::Abs(x))**2+([1]/2.)**2)*[2]*[3]**2", min(AxialPositions), max(AxialPositions))
        fSensCorrEW.FixParameter(0,AxialExt)
        fSensCorrEW.FixParameter(1,RingDiam)
        fSensCorrEW.FixParameter(2,PF*Eitotal**2)
        grASensGATE511keVEW.Fit(fSensCorrEW, "R")
        fSensCorrEW = grASensGATE511keVEW.GetFunction("fSensCorrEW")
        fSensCorrEW.SetLineColor(Color[0])
        fSensCorrEW.SetLineWidth(3)
        PhiFac.append( fSensCorrEW.GetParameter(3) )
        PhiFacErr.append( fSensCorrEW.GetParError(3) )
        RedChi2.append( fSensCorrEW.GetChisquare()/fSensCorrEW.GetNDF() )

        PSensEWPhi.append(fSensCorrEW.Eval(0.))

        MgrAxialSensEW.Add(grASensEW,'L')
        MgrAxialSensEW.Add(grASensGATE22NaEW,'P')
        MgrAxialSensEW.Add(grASensGATE511keVEW,'P')
        MgrAxialSensEW.Draw('A') 
        legendEW = ROOT.TLegend(0.25, 0.65, 0.99, 0.92)
        legendEW.AddEntry(grASensEW,'','L')
        legendEW.AddEntry(fSensCorrEW,'Correlated (fitting parameter #Phi)','L')
        legendEW.AddEntry(grASensGATE22NaEW,'','P')
        legendEW.AddEntry(grASensGATE511keVEW,'','P')
        legendEW.Draw('same')
        c1.Update()
        c1.SaveAs(output_folder[0]+'/AnalyticalSensitivity'+Scintillator+'EW.pdf')

        MgrAxialSensnoEW = ROOT.TMultiGraph()
        Title = Scintillator + '(no EW) ;axial position (mm);sensitivity (%)'
        MgrAxialSensnoEW.SetTitle(Title)
        MgrAxialSensnoEW.SetMaximum(65.)
        grASensnoEW = ROOT.TGraph( n_pos, pos, SensCorrnoEW )
        grASensnoEW.SetLineColor( Color[0] )
        grASensnoEW.SetMarkerColor( Color[0] )
        grASensnoEW.SetTitle(legend_text[0])
        grASens22NanoEW = ROOT.TGraph( n_pos, pos, Sens22NanoEW )
        grASens22NanoEW.SetLineColor( Color[1] )
        grASens22NanoEW.SetMarkerColor( Color[1] )
        grASens22NanoEW.SetTitle(legend_text[1])
        GateData = array( 'f', GateS22NanoEW[j])
        err = np.subtract(GateData, Sens22NanoEW)
        MAPEnoEW22Na.append( mape(GateData,Sens22NanoEW) )
        MAEnoEW22Na.append( mae(GateData,Sens22NanoEW) )
        STDEnoEW22Na.append( np.std(err) )
        grASensGATE22NanoEW = ROOT.TGraph( n_pos, AxialPositions, GateData)
        grASensGATE22NanoEW.SetLineColor( Color[1] )
        grASensGATE22NanoEW.SetMarkerColor( Color[1] )
        grASensGATE22NanoEW.SetTitle(legend_text[2])
        GateData = array( 'f', GateS511keVnoEW[j])
        err = np.subtract(GateData, SensCorrnoEW)
        MAPEnoEW511keV.append( mape(GateData,SensCorrnoEW) )
        MAEnoEW511keV.append( mae(GateData,SensCorrnoEW) )
        STDEnoEW511keV.append( np.std(err) )
        grASensGATE511keVnoEW = ROOT.TGraph( n_pos, AxialPositions, GateData )
        grASensGATE511keVnoEW.SetLineColor( Color[0] )
        grASensGATE511keVnoEW.SetMarkerColor( Color[0] )
        grASensGATE511keVnoEW.SetTitle(legend_text[3])
        MgrAxialSensnoEW.Add(grASensnoEW,'L')
        MgrAxialSensnoEW.Add(grASens22NanoEW,'L')
        MgrAxialSensnoEW.Add(grASensGATE22NanoEW,'P')
        MgrAxialSensnoEW.Add(grASensGATE511keVnoEW,'P')
        MgrAxialSensnoEW.Draw('A')
        legendnoEW = ROOT.TLegend(0.25, 0.65, 0.99, 0.92)
        legendnoEW.AddEntry(grASensnoEW,'','L')
        legendnoEW.AddEntry(grASens22NanoEW,'','L')
        legendnoEW.AddEntry(grASensGATE22NanoEW,'','P')
        legendnoEW.AddEntry(grASensGATE511keVnoEW,'','P')
        legendnoEW.Draw('same')
        c1.Update()
        c1.SaveAs(output_folder[0]+'/AnalyticalSensitivity'+Scintillator+'noEW.pdf')
    
    print ("Scintillator PSensnoEW511keVRef PSensnoEW511keV PSensnoEW511keVGATE PSensnoEW22Na PSensnoEW22NaGATE")
    for j, Scintillator in enumerate(Scintillators[start:stop], start=start):
        print(Scintillator+' '+str(round(PSensnoEW511keVRef[j],1))+" "+str(round(PSensnoEW511keV[j],1))+" "+str(round(PSensnoEW511keVGATE[j],1))+' '+str(round(PSensnoEW22Na[j],1))+" "+str(round(PSensnoEW22NaGATE[j],1)))
    print ("Scintillator Phifac RedChi2 PSensEWPE PSensEWPhi PSensEW511keVGATE")
    for j, Scintillator in enumerate(Scintillators[start:stop], start=start):
        print(Scintillator+' '+str(round(PhiFac[j],3))+u"\u00B1"+str(round(PhiFacErr[j],3))+" "+str(round(RedChi2[j],2))+" "+str(round(PSensEWPE[j],1))+" "+str(round(PSensEWPhi[j],1))+" "+str(round(PSensEW511keVGATE[j],1)))
    print ("Scintillator MAEnoEW22Na MAPEnoEW22Na STDEnoEW22Na MAEnoEW511keV MPAEnoEW511keV STDEnoEW511keV")
    for j, Scintillator in enumerate(Scintillators[start:stop], start=start):
        print(Scintillator+' '+str(round(MAEnoEW22Na[j],2))+" "+str(round(MAPEnoEW22Na[j],2))+" "+str(round(STDEnoEW22Na[j],2))+'  '+str(round(MAEnoEW511keV[j],2))+" "+str(round(MAPEnoEW511keV[j],2))+" "+str(round(STDEnoEW511keV[j],2)))
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

# --------------------------------------------------------------------------
if __name__ == '__main__':
    analysis_all_click()

