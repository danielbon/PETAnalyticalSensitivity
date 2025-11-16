#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 
# ./runSensitivityvsThicknessLYSOEW.py ./output/Sensitivity > SensitivityxThicknessLYSOEW.txt
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
    ThicknessLYSO_BaF2  = array( 'f', [ 25., 30., 35.] )
    eThicknessLYSO_BaF2 = array( 'f', [ 0., 0., 0.] )
    PeakSensLYSO_BaF2noEW  = array( 'f', [34.77, 41.50, 47.48] )
    MeanSensLYSO_BaF2noEW  = array( 'f', [ 23.85, 28.53, 32.70] )
    PeakSensLYSO_BaF2EW  = array( 'f', [ 6.84, 8.47, 9.90] )
    MeanSensLYSO_BaF2EW  = array( 'f', [ 4.23, 5.26, 6.19] )

    ThicknessLYSO_EJ232  = array( 'f', [ 25., 30., 35.] )
    eThicknessLYSO_EJ232 = array( 'f', [ 0., 0., 0.] )
    PeakSensLYSO_EJ232noEW  = array( 'f', [32.35, 38.77, 44.45] )
    MeanSensLYSO_EJ232noEW  = array( 'f', [ 22.00, 26.42, 30.41] )
    PeakSensLYSO_EJ232EW  = array( 'f', [ 6.57, 8.13, 9.53] )
    MeanSensLYSO_EJ232EW  = array( 'f', [4.01, 5.00, 5.90] )

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
    nLYSO_BaF2 = len(ThicknessLYSO_BaF2)
    nLYSO_EJ232 = len(ThicknessLYSO_EJ232)
    grThicknessSensitivity = ROOT.TMultiGraph()
    grThicknessSensitivity.SetTitle(';metascintillator thickness (mm);sensitivity (%)')
    grThicknessSensitivity.SetMaximum(12.)

    grPeakSensLYSO_BaF2EW = ROOT.TGraphErrors( nLYSO_BaF2, ThicknessLYSO_BaF2, PeakSensLYSO_BaF2EW, eThicknessLYSO_BaF2, eSens)
    grPeakSensLYSO_BaF2EW.SetTitle('Maximum sensitivity')
    grPeakSensLYSO_BaF2EW.SetMarkerColor( 12 )
    grPeakSensLYSO_BaF2EW.SetLineColor(12)
    fPeakSens = ROOT.TF1("fPeakSens", "[0] + x*[1]", min(ThicknessLYSO_BaF2), max(ThicknessLYSO_BaF2))
    grPeakSensLYSO_BaF2EW.Fit(fPeakSens)
    fPeakSens = grPeakSensLYSO_BaF2EW.GetFunction("fPeakSens")
    fPeakSens.SetLineColor(12)
    print ("Fit results: Linear Coefficient=",fPeakSens.GetParameter(0)," +/- ",fPeakSens.GetParError(0))
    print ("Fit results: Angular Coefficient=",fPeakSens.GetParameter(1)," +/- ",fPeakSens.GetParError(1))
    EqBGOThick[0] = (PeakSensBGOEW[0]-fPeakSens.GetParameter(0))/fPeakSens.GetParameter(1)
    print ("EqBGOThick for PeakSensLYSO_BaF2EW=",EqBGOThick)
    grPeakSensRefLYSO_BaF2EW = ROOT.TGraphErrors( 1, EqBGOThick, PeakSensBGOEW, eThicknessBGO, eSens)
    grPeakSensRefLYSO_BaF2EW.SetTitle('')
    grPeakSensRefLYSO_BaF2EW.SetMarkerColor( 1 )
    grPeakSensRefLYSO_BaF2EW.SetLineWidth( 0 )
    grPeakSensRefLYSO_BaF2EW.SetMarkerStyle(RefMarkerStyle)
    grPeakSensRefLYSO_BaF2EW.SetMarkerSize(RefMarkerSize)

    grMeanSensLYSO_BaF2EW = ROOT.TGraphErrors( nLYSO_BaF2, ThicknessLYSO_BaF2, MeanSensLYSO_BaF2EW, eThicknessLYSO_BaF2, eSens)
    grMeanSensLYSO_BaF2EW.SetTitle('Mean sensitivity')
    grMeanSensLYSO_BaF2EW.SetMarkerColor( 2 )
    grMeanSensLYSO_BaF2EW.SetLineColor(2)
    fMeanSens = ROOT.TF1("fMeanSens", "[0] + x*[1]", min(ThicknessLYSO_BaF2), max(ThicknessLYSO_BaF2))
    grMeanSensLYSO_BaF2EW.Fit(fMeanSens)
    fMeanSens = grMeanSensLYSO_BaF2EW.GetFunction("fMeanSens")
    fMeanSens.SetLineColor(2)
    print ("Fit results: Linear Coefficient=",fMeanSens.GetParameter(0)," +/- ",fMeanSens.GetParError(0))
    print ("Fit results: Angular Coefficient=",fMeanSens.GetParameter(1)," +/- ",fMeanSens.GetParError(1))
    EqBGOThick[0] = (MeanSensBGOEW[0]-fMeanSens.GetParameter(0))/fMeanSens.GetParameter(1)
    print ("EqBGOThick for MeanSensLYSO_BaF2EW=",EqBGOThick)
    grMeanSensRefLYSO_BaF2EW = ROOT.TGraphErrors( 1, EqBGOThick, MeanSensBGOEW, eThicknessBGO, eSens)
    grMeanSensRefLYSO_BaF2EW.SetTitle('')
    grMeanSensRefLYSO_BaF2EW.SetMarkerColor( 1 )
    grMeanSensRefLYSO_BaF2EW.SetLineWidth(0)
    grMeanSensRefLYSO_BaF2EW.SetMarkerStyle(RefMarkerStyle)
    grMeanSensRefLYSO_BaF2EW.SetMarkerSize(RefMarkerSize)

    grPeakSensLYSO_EJ232EW = ROOT.TGraphErrors( nLYSO_EJ232, ThicknessLYSO_EJ232, PeakSensLYSO_EJ232EW, eThicknessLYSO_EJ232, eSens)
    grPeakSensLYSO_EJ232EW.SetTitle('Maximum sensitivity')
    grPeakSensLYSO_EJ232EW.SetMarkerColor( 3 )
    grPeakSensLYSO_EJ232EW.SetLineColor(3)
    fPeakSens = ROOT.TF1("fPeakSens", "[0] + x*[1]", min(ThicknessLYSO_EJ232), max(ThicknessLYSO_EJ232))
    grPeakSensLYSO_EJ232EW.Fit(fPeakSens)
    fPeakSens = grPeakSensLYSO_EJ232EW.GetFunction("fPeakSens")
    fPeakSens.SetLineColor(3)
    print ("Fit results: Linear Coefficient=",fPeakSens.GetParameter(0)," +/- ",fPeakSens.GetParError(0))
    print ("Fit results: Angular Coefficient=",fPeakSens.GetParameter(1)," +/- ",fPeakSens.GetParError(1))
    EqBGOThick[0] = (PeakSensBGOEW[0]-fPeakSens.GetParameter(0))/fPeakSens.GetParameter(1)
    print ("EqBGOThick for PeakSensLYSO_EJ232EW=",EqBGOThick)
    grPeakSensRefLYSO_EJ232EW = ROOT.TGraphErrors( 1, EqBGOThick, PeakSensBGOEW, eThicknessBGO, eSens)
    grPeakSensRefLYSO_EJ232EW.SetTitle('')
    grPeakSensRefLYSO_EJ232EW.SetMarkerColor( 1 )
    grPeakSensRefLYSO_EJ232EW.SetLineWidth(0)
    grPeakSensRefLYSO_EJ232EW.SetMarkerStyle(RefMarkerStyle)
    grPeakSensRefLYSO_EJ232EW.SetMarkerSize(RefMarkerSize)

    grMeanSensLYSO_EJ232EW = ROOT.TGraphErrors( nLYSO_EJ232, ThicknessLYSO_EJ232, MeanSensLYSO_EJ232EW, eThicknessLYSO_EJ232, eSens)
    grMeanSensLYSO_EJ232EW.SetTitle('Mean sensitivity')
    grMeanSensLYSO_EJ232EW.SetMarkerColor( 4 )
    grMeanSensLYSO_EJ232EW.SetLineColor(4)
    fMeanSens = ROOT.TF1("fMeanSens", "[0] + x*[1]", min(ThicknessLYSO_EJ232), max(ThicknessLYSO_EJ232))
    grMeanSensLYSO_EJ232EW.Fit(fMeanSens)
    fMeanSens = grMeanSensLYSO_EJ232EW.GetFunction("fMeanSens")
    fMeanSens.SetLineColor(4)
    print ("Fit results: Linear Coefficient=",fMeanSens.GetParameter(0)," +/- ",fMeanSens.GetParError(0))
    print ("Fit results: Angular Coefficient=",fMeanSens.GetParameter(1)," +/- ",fMeanSens.GetParError(1))
    EqBGOThick[0] = (MeanSensBGOEW[0]-fMeanSens.GetParameter(0))/fMeanSens.GetParameter(1)
    print ("EqBGOThick for MeanSensLYSO_EJ232EW=",EqBGOThick)
    grMeanSensRefLYSO_EJ232EW = ROOT.TGraphErrors( 1, EqBGOThick, MeanSensBGOEW, eThicknessBGO, eSens)
    grMeanSensRefLYSO_EJ232EW.SetTitle('')
    grMeanSensRefLYSO_EJ232EW.SetMarkerColor( 1 )
    grMeanSensRefLYSO_EJ232EW.SetLineWidth(0)
    grMeanSensRefLYSO_EJ232EW.SetMarkerStyle(RefMarkerStyle)
    grMeanSensRefLYSO_EJ232EW.SetMarkerSize(RefMarkerSize)

    grThicknessSensitivity.Add(grPeakSensLYSO_BaF2EW)
    grThicknessSensitivity.Add(grMeanSensLYSO_BaF2EW)
    grThicknessSensitivity.Add(grPeakSensRefLYSO_BaF2EW)
    grThicknessSensitivity.Add(grMeanSensRefLYSO_BaF2EW)
    grThicknessSensitivity.Add(grPeakSensLYSO_EJ232EW)
    grThicknessSensitivity.Add(grMeanSensLYSO_EJ232EW)
    grThicknessSensitivity.Add(grPeakSensRefLYSO_EJ232EW)
    grThicknessSensitivity.Add(grMeanSensRefLYSO_EJ232EW)
    grThicknessSensitivity.Draw( 'AP' )

    legend = ROOT.TLegend(0.17, 0.64, 0.67, 0.94)
    legend.AddEntry(grPeakSensLYSO_BaF2EW ,'Peak - LYSO/BaF2')
    legend.AddEntry(grMeanSensLYSO_BaF2EW ,'Mean - LYSO/BaF2')
    legend.AddEntry(grPeakSensLYSO_EJ232EW ,'Peak - LYSO/EJ232')
    legend.AddEntry(grMeanSensLYSO_EJ232EW ,'Mean - LYSO/EJ232')
    legend.AddEntry(grPeakSensRefLYSO_BaF2EW ,'Reference (BGO)')
    legend.Draw('same')
    c1.Update() 
    c1.SaveAs("SensitivityxThicknessLYSOEW.pdf")
    c1.SaveAs("SensitivityxThicknessLYSOEW.png")
    r = True
    return r

# --------------------------------------------------------------------------
if __name__ == '__main__':
    analysis_all_click()

