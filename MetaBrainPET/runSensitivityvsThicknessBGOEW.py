#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 
# ./runSensitivityvsThicknessBGOEW.py ./output/Sensitivity > SensitivityxThicknessBGOEWnew.txt
from __future__ import print_function

import uproot
import click
import gatetools as gt
import os
import numpy as np

from box import Box
from statistics import mean

import ROOT
#from ROOT import TCanvas, TGraph
#from ROOT import gROOT
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
    c1.SetBottomMargin(0.14)
    c1.SetLeftMargin(0.10)
    
    RefMarkerStyle = 43
    RefMarkerSize = 5
    ThicknessBGO_BaF2  = array( 'f', [ 15., 20., 25.] )
    eThicknessBGO_BaF2 = array( 'f', [ 0., 0., 0.] )
    PeakSensBGO_BaF2noEW  = array( 'f', [ 21.03, 29.19, 36.63 ] )
    MeanSensBGO_BaF2noEW  = array( 'f', [ 14.29, 19.84, 24.95 ] )
    PeakSensBGO_BaF2EW  = array( 'f', [ 4.78, 7.11, 9.29 ] )
    MeanSensBGO_BaF2EW  = array( 'f', [ 2.86, 4.30, 5.67 ] )

    ThicknessBGO_EJ232  = array( 'f', [ 15., 20., 25.] )
    eThicknessBGO_EJ232 = array( 'f', [ 0., 0., 0.] )
    PeakSensBGO_EJ232noEW  = array( 'f', [20.38, 28.28, 35.44] )
    MeanSensBGO_EJ232noEW  = array( 'f', [ 12.34, 17.13, 21.53] )
    PeakSensBGO_EJ232EW  = array( 'f', [ 5.11, 7.53, 9.78] )
    MeanSensBGO_EJ232EW  = array( 'f', [3.01, 4.49, 5.88] )

    PeakSensLYSOnoEW  = array( 'f', [33.72] )
    MeanSensLYSOnoEW  = array( 'f', [22.96] )
    PeakSensLYSOEW  = array( 'f', [8.47] )
    MeanSensLYSOEW  = array( 'f', [5.14] )
    ThicknessLYSO  = array( 'f', [20.] )
    eThicknessLYSO = array( 'f', [0.] )
    EqLYSOThick = array( 'f', [0.] )

    PeakSensBGOnoEW  = array( 'f', [27.83] )
    MeanSensBGOnoEW  = array( 'f', [18.71] )
    PeakSensBGOEW  = array( 'f', [8.76] )
    MeanSensBGOEW  = array( 'f', [5.18] )
    ThicknessBGO  = array( 'f', [15.] )
    eThicknessBGO = array( 'f', [0.] )
    EqBGOThick = array( 'f', [0.] )

    eSens = array( 'f', [ 0.01, 0.01, 0.01, 0.01 ] ) 
    nBGO_BaF2 = len(ThicknessBGO_BaF2)
    nBGO_EJ232 = len(ThicknessBGO_EJ232)
    grThicknessSensitivity = ROOT.TMultiGraph()
    grThicknessSensitivity.SetTitle(';metascintillator thickness (mm); ')
    grThicknessSensitivity.SetMaximum(12.)

    grPeakSensBGO_BaF2EW = ROOT.TGraphErrors( nBGO_BaF2, ThicknessBGO_BaF2, PeakSensBGO_BaF2EW, eThicknessBGO_BaF2, eSens)
    grPeakSensBGO_BaF2EW.SetTitle('Maximum sensitivity')
    grPeakSensBGO_BaF2EW.SetMarkerColor( 12 )
    grPeakSensBGO_BaF2EW.SetLineColor(12)
    fPeakSens = ROOT.TF1("fPeakSens", "[0] + x*[1]", min(ThicknessBGO_BaF2), max(ThicknessBGO_BaF2))
    grPeakSensBGO_BaF2EW.Fit(fPeakSens)
    fPeakSens = grPeakSensBGO_BaF2EW.GetFunction("fPeakSens")
    fPeakSens.SetLineColor(12)
    print ("Fit results: Linear Coefficient=",fPeakSens.GetParameter(0)," +/- ",fPeakSens.GetParError(0))
    print ("Fit results: Angular Coefficient=",fPeakSens.GetParameter(1)," +/- ",fPeakSens.GetParError(1))
    EqBGOThick[0] = (PeakSensBGOEW[0]-fPeakSens.GetParameter(0))/fPeakSens.GetParameter(1)
    print ("EqBGOThick for PeakSensBGO_BaF2EW=",EqBGOThick)
    grPeakSensRefBGO_BaF2EW = ROOT.TGraphErrors( 1, EqBGOThick, PeakSensBGOEW, eThicknessBGO, eSens)
    grPeakSensRefBGO_BaF2EW.SetTitle('')
    grPeakSensRefBGO_BaF2EW.SetMarkerColor(1)
    grPeakSensRefBGO_BaF2EW.SetLineWidth(0)
    grPeakSensRefBGO_BaF2EW.SetMarkerStyle(RefMarkerStyle)
    grPeakSensRefBGO_BaF2EW.SetMarkerSize(RefMarkerSize)

    grMeanSensBGO_BaF2EW = ROOT.TGraphErrors( nBGO_BaF2, ThicknessBGO_BaF2, MeanSensBGO_BaF2EW, eThicknessBGO_BaF2, eSens)
    grMeanSensBGO_BaF2EW.SetTitle('Mean sensitivity')
    grMeanSensBGO_BaF2EW.SetMarkerColor( 2 )
    grMeanSensBGO_BaF2EW.SetLineColor(2)
    fMeanSens = ROOT.TF1("fMeanSens", "[0] + x*[1]", min(ThicknessBGO_BaF2), max(ThicknessBGO_BaF2))
    grMeanSensBGO_BaF2EW.Fit(fMeanSens)
    fMeanSens = grMeanSensBGO_BaF2EW.GetFunction("fMeanSens")
    fMeanSens.SetLineColor(2)
    print ("Fit results: Linear Coefficient=",fMeanSens.GetParameter(0)," +/- ",fMeanSens.GetParError(0))
    print ("Fit results: Angular Coefficient=",fMeanSens.GetParameter(1)," +/- ",fMeanSens.GetParError(1))
    EqBGOThick[0] = (MeanSensBGOEW[0]-fMeanSens.GetParameter(0))/fMeanSens.GetParameter(1)
    print ("EqBGOThick for MeanSensBGO_BaF2EW=",EqBGOThick)
    grMeanSensRefBGO_BaF2EW = ROOT.TGraphErrors( 1, EqBGOThick, MeanSensBGOEW, eThicknessBGO, eSens)
    grMeanSensRefBGO_BaF2EW.SetTitle('')
    grMeanSensRefBGO_BaF2EW.SetMarkerColor(1)
    grMeanSensRefBGO_BaF2EW.SetLineWidth(0)
    grMeanSensRefBGO_BaF2EW.SetMarkerStyle(RefMarkerStyle)
    grMeanSensRefBGO_BaF2EW.SetMarkerSize(RefMarkerSize)

    grPeakSensBGO_EJ232EW = ROOT.TGraphErrors( nBGO_EJ232, ThicknessBGO_EJ232, PeakSensBGO_EJ232EW, eThicknessBGO_EJ232, eSens)
    grPeakSensBGO_EJ232EW.SetTitle('Maximum sensitivity')
    grPeakSensBGO_EJ232EW.SetMarkerColor( 3 )
    grPeakSensBGO_EJ232EW.SetLineColor(3)
    fPeakSens = ROOT.TF1("fPeakSens", "[0] + x*[1]", min(ThicknessBGO_EJ232), max(ThicknessBGO_EJ232))
    grPeakSensBGO_EJ232EW.Fit(fPeakSens)
    fPeakSens = grPeakSensBGO_EJ232EW.GetFunction("fPeakSens")
    fPeakSens.SetLineColor(3)
    print ("Fit results: Linear Coefficient=",fPeakSens.GetParameter(0)," +/- ",fPeakSens.GetParError(0))
    print ("Fit results: Angular Coefficient=",fPeakSens.GetParameter(1)," +/- ",fPeakSens.GetParError(1))
    EqBGOThick[0] = (PeakSensBGOEW[0]-fPeakSens.GetParameter(0))/fPeakSens.GetParameter(1)
    print ("EqBGOThick for PeakSensBGO_EJ232EW=",EqBGOThick)
    grPeakSensRefBGO_EJ232EW = ROOT.TGraphErrors( 1, EqBGOThick, PeakSensBGOEW, eThicknessBGO, eSens)
    grPeakSensRefBGO_EJ232EW.SetTitle('')
    grPeakSensRefBGO_EJ232EW.SetMarkerColor(1)
    grPeakSensRefBGO_EJ232EW.SetLineWidth(0)
    grPeakSensRefBGO_EJ232EW.SetMarkerStyle(RefMarkerStyle)
    grPeakSensRefBGO_EJ232EW.SetMarkerSize(RefMarkerSize)

    grMeanSensBGO_EJ232EW = ROOT.TGraphErrors( nBGO_EJ232, ThicknessBGO_EJ232, MeanSensBGO_EJ232EW, eThicknessBGO_EJ232, eSens)
    grMeanSensBGO_EJ232EW.SetTitle('Mean sensitivity')
    grMeanSensBGO_EJ232EW.SetMarkerColor( 4 )
    grMeanSensBGO_EJ232EW.SetLineColor(4)
    fMeanSens = ROOT.TF1("fMeanSens", "[0] + x*[1]", min(ThicknessBGO_EJ232), max(ThicknessBGO_EJ232))
    grMeanSensBGO_EJ232EW.Fit(fMeanSens)
    fMeanSens = grMeanSensBGO_EJ232EW.GetFunction("fMeanSens")
    fMeanSens.SetLineColor(4)
    print ("Fit results: Linear Coefficient=",fMeanSens.GetParameter(0)," +/- ",fMeanSens.GetParError(0))
    print ("Fit results: Angular Coefficient=",fMeanSens.GetParameter(1)," +/- ",fMeanSens.GetParError(1))
    EqBGOThick[0] = (MeanSensBGOEW[0]-fMeanSens.GetParameter(0))/fMeanSens.GetParameter(1)
    print ("EqBGOThick for MeanSensBGO_EJ232EW=",EqBGOThick)
    grMeanSensRefBGO_EJ232EW = ROOT.TGraphErrors( 1, EqBGOThick, MeanSensBGOEW, eThicknessBGO, eSens)
    grMeanSensRefBGO_EJ232EW.SetTitle('')
    grMeanSensRefBGO_EJ232EW.SetMarkerColor(1)
    grMeanSensRefBGO_EJ232EW.SetLineWidth(0)
    grMeanSensRefBGO_EJ232EW.SetMarkerStyle(RefMarkerStyle)
    grMeanSensRefBGO_EJ232EW.SetMarkerSize(RefMarkerSize)

    grThicknessSensitivity.Add(grPeakSensBGO_BaF2EW)
    grThicknessSensitivity.Add(grMeanSensBGO_BaF2EW)
    grThicknessSensitivity.Add(grPeakSensRefBGO_BaF2EW)
    grThicknessSensitivity.Add(grMeanSensRefBGO_BaF2EW)
    grThicknessSensitivity.Add(grPeakSensBGO_EJ232EW)
    grThicknessSensitivity.Add(grMeanSensBGO_EJ232EW)
    grThicknessSensitivity.Add(grPeakSensRefBGO_EJ232EW)
    grThicknessSensitivity.Add(grMeanSensRefBGO_EJ232EW)
    grThicknessSensitivity.Draw( 'AP' )

    legend = ROOT.TLegend(0.17, 0.64, 0.67, 0.94)
    legend.AddEntry(grPeakSensBGO_BaF2EW ,'Peak - BGO/BaF2')
    legend.AddEntry(grMeanSensBGO_BaF2EW ,'Mean - BGO/BaF2')
    legend.AddEntry(grPeakSensBGO_EJ232EW ,'Peak - BGO/EJ232')
    legend.AddEntry(grMeanSensBGO_EJ232EW ,'Mean - BGO/EJ232')
    legend.AddEntry(grPeakSensRefBGO_BaF2EW ,'Reference (BGO)')
    legend.Draw('same')
    c1.Update() 
    c1.SaveAs("SensitivityxThicknessBGOEW.pdf")
    c1.SaveAs("SensitivityxThicknessBGOEW.png")
    r = True
    return r

# --------------------------------------------------------------------------
if __name__ == '__main__':
    analysis_all_click()

