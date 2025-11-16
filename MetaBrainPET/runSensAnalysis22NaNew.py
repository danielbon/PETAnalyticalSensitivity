#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 
# ../../../runSensAnalysis22Na.py LYSO-Proteus/ BGO/ BGO_EJ232_15/ BGO_EJ232_20/ BGO_EJ232_25/ LYSO-Proteus_EJ232_25/ LYSO-Proteus_EJ232_30/ LYSO-Proteus_EJ232_35/ > ./output22NaLYSO-BGO_EJ232.txt 
# ../../../runSensAnalysis22Na.py LYSO-Proteus/ BGO/ BGO_BaF2_15/ BGO_BaF2_20/ BGO_BaF2_25/ LYSO-Proteus_BaF2_25/ LYSO-Proteus_BaF2_30/ LYSO-Proteus_BaF2_35/ > ./output22NaLYSO-BGO_BaF2.txt
# ../../../runSensAnalysis22Na.py LYSO-Proteus/ BGO/ BGO_EJ232_20/ BGO_EJ232_25/ LYSO-Proteus_EJ232_25/ LYSO-Proteus_EJ232_30/ LYSO-Proteus_EJ232_35/ > ./output22NaLYSO-BGO_EJ232.txt 
# ../../../runSensAnalysis22Na.py LYSO-Proteus/ BGO/ BGO_BaF2_20/ BGO_BaF2_25/ LYSO-Proteus_BaF2_25/ LYSO-Proteus_BaF2_30/ LYSO-Proteus_BaF2_35/ > ./output22NaLYSO-BGO_BaF2.txt
# ../../../runSensAnalysis22NaNew.py LYSO-Proteus/ BGO/ BGO_BaF2_15/ BGO_BaF2_20/ BGO_BaF2_25/ LYSO-Proteus_BaF2_25/ LYSO-Proteus_BaF2_30/ LYSO-Proteus_BaF2_35/ BGO_EJ232_15/ BGO_EJ232_20/ BGO_EJ232_25/ LYSO-Proteus_EJ232_25/ LYSO-Proteus_EJ232_30/ LYSO-Proteus_EJ232_35/ > ./output22NaAllnew.txt
# ../../../runSensAnalysis511keV.py LYSO-Proteus/ BGO/ BGO_EJ232_25/ LYSO-Proteus_EJ232_30/ > ./output511keVLYSO-BGO_EJ232.txt 
# ../../../runSensAnalysis511keV.py LYSO-Proteus/ BGO/ BGO_BaF2_25/ LYSO-Proteus_BaF2_30/ > ./output511keVLYSO-BGO_BaF2.txt
#
#
from __future__ import print_function

import uproot
import click
import gatetools as gt
import os
import numpy as np
from box import Box
from statistics import mean
import ROOT
from math import sin
from array import array

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('output_folders', nargs=-1, required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=True))
@gt.add_options(gt.common_options)
def analysis_all_click(output_folders, **kwargs):
    analyse_all_folders(output_folders)

def analyse_all_folders(output_folders, **kwargs):
    FontSize = 0.034
    ROOT.gStyle.SetLabelSize(FontSize,"xyz")
    ROOT.gStyle.SetTitleSize(FontSize,"xyz")
    ROOT.gStyle.SetTextSize(FontSize)
    ROOT.gStyle.SetLegendTextSize(FontSize*6)
    ROOT.gStyle.SetLineWidth(2)
    ROOT.gStyle.SetLineStyle(1)
    ROOT.gStyle.SetMarkerSize(2)
    ROOT.gStyle.SetMarkerStyle(20)
    ROOT.gROOT.ForceStyle()
    ########### USER SECTION - These parameters must be verified before running the script#####
    BF = 0.9 # Branching fraction
    # Total activity times the Branching fraction (for sodium-22 is 0.9. Fluorine-18 is 0.97)
    a = 370000*BF  # in Bq
    print(f'Total activity             = {a / 1e6:.2f} MBq')
    n_folders = len(output_folders)
    pos_folders = ['./-75/', './-60/', './-45/', './-30/', './-15/', './0/', './15/', './30/', './45/', './60/', './75/']
    #legend_text_all = ['LYSO-20 mm', 'BGO-15 mm', 'BGO/'+'BaF2'+'-20 mm',  'BGO/'+'BaF2'+'-25 mm',  'LYSO/'+'BaF2'+'-25 mm',  'LYSO/'+'BaF2'+'-30 mm', 'LYSO/'+'BaF2'+'-35 mm', 'BGO/'+'EJ232'+'-20 mm',  'BGO/'+'EJ232'+'-25 mm',  'LYSO/'+'EJ232'+'-25 mm',  'LYSO/'+'EJ232'+'-30 mm', 'LYSO/'+'EJ232'+'-35 mm']
    legend_text_all = ['LYSO-20 mm', 'BGO-15 mm', 'BGO/'+'BaF2'+'-15 mm', 'BGO/'+'BaF2'+'-20 mm',  'BGO/'+'BaF2'+'-25 mm',  'LYSO/'+'BaF2'+'-25 mm',  'LYSO/'+'BaF2'+'-30 mm', 'LYSO/'+'BaF2'+'-35 mm', 'BGO/'+'EJ232'+'-15 mm','BGO/'+'EJ232'+'-20 mm',  'BGO/'+'EJ232'+'-25 mm',  'LYSO/'+'EJ232'+'-25 mm',  'LYSO/'+'EJ232'+'-30 mm', 'LYSO/'+'EJ232'+'-35 mm']
    ########### END OF USER SECTION ###########################################################
    n_pos = len(pos_folders)
    #LineColor=[1,2,3,3,4,4,4,9,9,28,28,28]
    #LineStyle=[2,2,3,7,3,7,9,3,7,3,7,9]
    LineColor=[1,2,3,3,3,4,4,4,9,9,9,28,28,28]
    LineStyle=[2,2,3,7,9,3,7,9,3,7,9,3,7,9]
#    LineColor=[1,2,3,4,6,7,8,9,12,15,28,40,38]
#    LineStyle=[1,2,3,4,6,7,8,9,12,15,28,40,38]
    #MgrAxialSensAll = ROOT.TMultiGraph()
    #MgrAxialSensAll.SetTitle(';axial position (mm);sensitivity (%)')
    #MgrAxialSensAll.SetMaximum(85.) 
    MgrAxialSensEW = ROOT.TMultiGraph()
    MgrAxialSensEW.SetTitle(';axial position (mm);')
    #MgrAxialSensEW.SetMinimum(0.)
    MgrAxialSensEW.SetMaximum(11.)
    #MgrAxialSensEW.GetXaxis().SetLimits(-85.,90.)
    MgrAxialSensnoEW = ROOT.TMultiGraph()
    MgrAxialSensnoEW.SetTitle(';axial position (mm);sensitivity (%)')
    #MgrAxialSensnoEW.SetMinimum(0.)
    MgrAxialSensnoEW.SetMaximum(50.)
    #MgrAxialSensnoEW.GetXaxis().SetLimits(-85.,90.)
    grASensEW = []
    grASensnoEW = []
    MaxSensnoEW, MeanSensnoEW, MaxSensEW, MeanSensEW = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )
    # loop over folders
    for index, output_folder in enumerate(output_folders, start=0):
    #for output_folder in output_folders:
        print(f' ')
        print(f'########## Crystal material and thickness (in mm) = {output_folder[0:-1]} #########')
        pos, SensitivitynoEW = array( 'd' ), array( 'd' )
        posEW, SensitivityEW = array( 'd' ), array( 'd' )
        for pos_folder in pos_folders:
            print(f'########## Axial position in mm (or pos_folder) = {pos_folder[2:-1]} #########')
            fn = os.path.join(output_folder+pos_folder, 'output.root')
            f = uproot.open(fn)
            print(f'########## Total Counts (no energy window) #########')
            data = analysis(f,'Coincidences', 'Delay')
            Rt = data.trues_count_rate
            Rtot = data.prompts_count_rate
            Rsc = data.scatter_count_rate
            necr = Rt ** 2 / Rtot
            sf = Rsc / (Rt + Rsc) * 100
            Sensitivity = Rt / a * 100
            print(f'Simulation duration        = {data.duration} sec')
            print(f'Sensitivity                = {Sensitivity:.2f} %')
            print(f'ScatterFraction            = {sf:.2f} %   ')
            print(f'NECR                       = {necr:.2f} cps')
            pos.append( float(pos_folder[2:-1]) )
            SensitivitynoEW.append( Sensitivity )
            print(f'########## Counts using Energy Window (EW) from 350 keV to 650 keV #########')
            data = analysis(f,'Coincidences300', 'Delay300')
            Rt = data.trues_count_rate
            Rtot = data.prompts_count_rate
            Rsc = data.scatter_count_rate
            necr = Rt ** 2 / Rtot
            sf = Rsc / (Rt + Rsc) * 100
            Sensitivity = Rt / a * 100
            print(f'Simulation duration        = {data.duration} sec')
            print(f'Sensitivity                = {Sensitivity:.2f} %')
            print(f'ScatterFraction            = {sf:.2f} %   ')
            print(f'NECR                       = {necr:.2f} cps')
            posEW.append( float(pos_folder[2:-1]) )
            SensitivityEW.append( Sensitivity )
        pos, SensitivitynoEW = zip(*sorted(zip(pos, SensitivitynoEW)))
        pos = np.asarray(pos)
        SensitivitynoEW = np.asarray(SensitivitynoEW)
        posEW, SensitivityEW = zip(*sorted(zip(posEW, SensitivityEW)))
        posEW = np.asarray(posEW)
        SensitivityEW = np.asarray(SensitivityEW)
        print(f'Position (mm)  Sensitivity without energy window') 
        for i in range( n_pos ):
           print(' %.1f   %.2f ' % (pos[i],SensitivitynoEW[i]))
        print(f'Position (mm)  Sensitivity with energy window (350 - 650 keV)') 
        for i in range( n_pos ):
           print(' %.1f   %.2f ' % (posEW[i],SensitivityEW[i]))
        MaxSensnoEW.append(max(SensitivitynoEW))
        MeanSensnoEW.append(mean(SensitivitynoEW))
        MaxSensEW.append(max(SensitivityEW))
        MeanSensEW.append(mean(SensitivityEW))
        grASensEW.append(ROOT.TGraph( n_pos, posEW, SensitivityEW ))
        grASensEW[index].SetLineColor( LineColor[index] )
        grASensEW[index].SetLineStyle(LineStyle[index])
        grASensEW[index].SetMarkerColor( LineColor[index] )
        grASensEW[index].SetTitle(legend_text_all[index])
        grASensnoEW.append( ROOT.TGraph( n_pos, pos, SensitivitynoEW ))
        grASensnoEW[index].SetLineColor( LineColor[index] )
        grASensnoEW[index].SetLineStyle(LineStyle[index])
        grASensnoEW[index].SetMarkerColor( LineColor[index] )
        grASensnoEW[index].SetTitle(legend_text_all[index])
        MgrAxialSensnoEW.Add(grASensnoEW[index])
        MgrAxialSensEW.Add(grASensEW[index])

    print(f' ')
    print("Crystal")
    for i in range( n_folders ):
        print(f'{output_folders[i]},', end=" ")
    print("")
    print("MaxSensnoEW")    
    for i in range( n_folders ):
        print(f'{MaxSensnoEW[i]:.2f},', end=" ")
    print("")
    print("MeanSensnoEW")    
    for i in range( n_folders ):
        print(f'{MeanSensnoEW[i]:.2f},', end=" ")
    print("")
    print("MaxSensEW")    
    for i in range( n_folders ):
        print(f'{MaxSensEW[i]:.2f},', end=" ")
    print("")
    print("MeanSensEW")    
    for i in range( n_folders ):
        print(f'{MeanSensEW[i]:.2f},', end=" ")
    print("")

    c1 = ROOT.TCanvas( 'c1', '', 1500, 1200 )
    c1.SetGridx(1)
    c1.SetGridy(0)
    c1.SetRightMargin(0.02)
    c1.SetTopMargin(0.02)
    c1.SetBottomMargin(0.10)
    c1.SetLeftMargin(0.10)

    pad1 = ROOT.TPad("pad1","pad1" ,0.0 ,0.0 ,0.5 ,0.85)
    #pad1.SetBottomMargin(0)
    pad1.SetTopMargin(0.02)
    pad1.SetRightMargin(0.02)
    pad1.Draw()
    pad1.cd()
    MgrAxialSensnoEW.Draw('APL')
    latex1 = ROOT.TLatex()
    latex1.SetNDC()
    #latex.SetTextSize(0.06)
    latex1.DrawText(0.7 ,0.83 , "no EW")

    #legendnoEW = ROOT.TLegend(0.15, 0.65, 0.99, 0.99)
    #legendnoEW.SetNColumns(2)
    #for index, output_folder in enumerate(output_folders, start=0):
    #    legendnoEW.AddEntry(grASensnoEW[index] ,legend_text_all[index])
    #legendnoEW.Draw('same')
    #c1.Update()
    #c1.SaveAs(output_folder+"../sensitivitynoEW.pdf")
    c1.cd()
    pad2 = ROOT.TPad("pad2","pad2" ,0.5 ,0.0 ,1.0 ,0.85)
    #pad2.SetBottomMargin(0)
    pad2.SetLeftMargin(0.05)
    pad2.SetTopMargin(0.02)
    pad2.SetRightMargin(0.02)
    pad2.Draw()
    pad2.cd()
    MgrAxialSensEW.Draw('APL')
    latex2 = ROOT.TLatex()
    latex2.SetNDC()
    latex2.DrawText(0.7 ,0.83 , "EW")

    c1.cd()
    pad3 = ROOT.TPad("pad3","pad3" ,0.01 ,0.855 ,1.0 ,1.0)
    pad3.SetBottomMargin(0.02)
    #pad3.SetTopMargin(0)
    pad3.Draw()
    pad3.cd()
    legendEW = ROOT.TLegend(0., 0., 0.99, 0.99)
    legendEW.SetNColumns(4)
    for index, output_folder in enumerate(output_folders, start=0):
        legendEW.AddEntry(grASensEW[index] ,legend_text_all[index])
    legendEW.Draw('same')

    c1.Update()
    c1.SaveAs(output_folder+"../sensitivityAllnew.pdf")

    r = True
    return r

def get_pet_counts_plus(root_file, coincidences_evts, delay_evts):
    data = Box()
    coincidences = root_file[coincidences_evts]
    data.prompts_count = coincidences.num_entries
    delays = root_file[delay_evts]
    data.delays_count = delays.num_entries
    # Root tree analysis principle:
    # 1) consider a "Tree" (such as root_file['Coincidences'])
    # 2) convert the branches to numpy array
    # 3) use mask to select only some part of the tree
    # Compute the number of Random Coincidences
    # (when the two singles in coincidence came from two different events or are noise)
    event_id1 = coincidences['eventID1'].array(library='numpy')
    event_id2 = coincidences['eventID2'].array(library='numpy')
    # Warning noise events have eventID == -2
    # https://opengate.readthedocs.io/en/latest/digitizer_and_detector_modeling.html?#noise
    mask_trues = (event_id1 == event_id2) & (event_id1 >= 0)
    data.randoms_count = data.prompts_count - len(event_id1[mask_trues])
    # Scattered coincidences (among the true)
    scatter_compton1 = coincidences['comptonPhantom1'].array(library='numpy')
    scatter_compton2 = coincidences['comptonPhantom2'].array(library='numpy')
    scatter_rayleigh1 = coincidences['RayleighPhantom1'].array(library='numpy')
    scatter_rayleigh2 = coincidences['RayleighPhantom2'].array(library='numpy')
    # if any of the value is greater than zero, it means a compton or a rayleigh occurs, so this is a scattered event
    mask_scatter = (scatter_compton1 > 0) | (scatter_compton2 > 0) | (scatter_rayleigh1 > 0) | (scatter_rayleigh2 > 0)
    data.scatter_count = len(scatter_compton1[mask_trues & mask_scatter])
    # Remaining true events
    data.trues_count = data.prompts_count - data.randoms_count - data.scatter_count
    return data

def get_pet_data(root_file):
    data = Box()
    try:
        data_pet = root_file['pet_data']
        data.total_nb_primaries = data_pet['total_nb_primaries'].array(library='numpy')[0]
        data.latest_event_ID = data_pet['latest_event_ID'].array(library='numpy')[0]
        data.stop_time_sec = data_pet['stop_time_sec'].array(library='numpy')[0]
        data.start_time_sec = data_pet['start_time_sec'].array(library='numpy')[0]
        return data
    except:
        return False

def analysis(root_file, coincidences_tree, delay_tree):
    # Compute the types of events
    data = get_pet_counts_plus(root_file, coincidences_tree, delay_tree)
    # Compute the rate -> divide the counts per seconds (if available)
    d = get_pet_data(root_file)
    if not d:
        print('Need the acquisition time to compute the NECR')
        exit(0)
    data.duration = d.stop_time_sec - d.start_time_sec
    print(f'Prompt events    {data.prompts_count} \t(all coincidences)')
    print(f'Delayed events   {data.delays_count}  \t(approximation of the number of random events)')
    print(f'Random events    {data.randoms_count} \t(including noise event)')
    print(f'Scattered events {data.scatter_count}')
    print(f'Trues events     {data.trues_count}')
    data.prompts_count_rate = data.prompts_count / data.duration
    data.delays_count_rate = data.delays_count / data.duration
    data.randoms_count_rate = data.randoms_count / data.duration
    data.scatter_count_rate = data.scatter_count / data.duration
    data.trues_count_rate = data.trues_count / data.duration
    print(f'Prompt rate     {data.prompts_count_rate:.0f} cps')
    print(f'Delayed rate    {data.delays_count_rate:.0f} cps')
    print(f'Random rate     {data.randoms_count_rate:.0f} cps')
    print(f'Scattered rate  {data.scatter_count_rate:.0f} cps')
    print(f'Trues rate      {data.trues_count_rate:.0f} cps')
    return data

# --------------------------------------------------------------------------
if __name__ == '__main__':
    analysis_all_click()
