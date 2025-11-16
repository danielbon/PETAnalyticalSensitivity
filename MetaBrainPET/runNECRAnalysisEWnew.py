#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 
# ../../runNECRAnalysisEW.py LYSO-Proteus/ BGO/ BGO_BaF2_25/ BGO_EJ232_25/ LYSO-Proteus_BaF2_25/ LYSO-Proteus_EJ232_25/ LYSO-Proteus_BaF2_30/ LYSO-Proteus_EJ232_30/ > ./nec-output_all.txt
# ../../runNECRAnalysisEW.py LYSO-Proteus/ BGO/ BGO_BaF2_23.2/ BGO_EJ232_22.5/ LYSO-Proteus_BaF2_29.7/ LYSO-Proteus_EJ232_31.1/ > ./nec-output_all.txt
# ../../runNECRAnalysisEWnew.py LYSO-Proteus/ BGO/ BGO_BaF2_23.7/ BGO_EJ232_22.7/ LYSO-Proteus_BaF2_31.1/ LYSO-Proteus_EJ232_32.3/ > ./nec-output_allnew.txt
# ../../runNECRAnalysisEW.py LYSO-Proteus/ BGO/ BGO_BaF2_23.2/ BGO_EJ232_22.5/ > ./nec-output_all.txt
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
    FontSize = 0.04
    ROOT.gStyle.SetLabelSize(FontSize,"xyz")
    ROOT.gStyle.SetTitleSize(FontSize,"xyz")
    ROOT.gStyle.SetTextSize(FontSize)
    ROOT.gStyle.SetLegendTextSize(6*FontSize)
    ROOT.gStyle.SetLineWidth( 1 )
    ROOT.gStyle.SetLineStyle( 1 )
    ROOT.gStyle.SetMarkerSize( 2 )
    ROOT.gStyle.SetMarkerStyle( 20 )
    ROOT.gROOT.ForceStyle()
    Color=[1,2,3,4,6,9,12,15,28]
    gt.logging_conf(**kwargs)
    n_folders = len(output_folders)
    activity_folders = ['./10/', './50/', './100/', './300/', './400/', './500/', './750/', './1000/', './1500/', './2000/', './2500/', './3000/']
    n_activities = len(activity_folders)
    maximum = 140000.
    maxToF = 10.5
    grNECR = ROOT.TMultiGraph()
    grNECR.SetTitle(';activity (MBq);NECR (cps)')
    #grNECR.SetMinimum(0.)
    grNECR.SetMaximum(maximum)
    grNECREW = []
    grNECRToF = ROOT.TMultiGraph()
    grNECRToF.SetTitle(';activity (MBq);')
    #grNECR.SetMinimum(0.)
    grNECRToF.SetMaximum(maximum*maxToF)
    grNECREWToF = []
    #grNECRnoEW = []
    MaxNECRnoEW, MeanNECRnoEW, MaxNECREW, MeanNECREW = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )
    source_vol = 0.16*0.16*70*3.1416
    for index, output_folder in enumerate(output_folders, start=0):
        print(f' ')
        print(f'########## Crystal material and thickness (in mm) = {output_folder[0:-1]} #########')
        #activity, NECRnoEW = array( 'd' ), array( 'd' )
        activityEW, NECREW = array( 'd' ), array( 'd' )
        NECREWToF = array( 'd' )
        for activity_folder in activity_folders:
            print(f'########## Activity in MBq (or activity_folder) = {activity_folder[2:-1]} #########')
            BF = 1.
            # Total activity = activity times the branching ratio
            a = float(activity_folder[2:-1])*BF*1e6  # in Bq
            print(f'Total activity = activity times the branching ratio = {a / 1e6:.2f} MBq')
            fn = os.path.join(output_folder+activity_folder, 'output.root')
            f = uproot.open(fn)
            #print(f'########## Total Counts (no energy window) #########')
            #data = analysis(f,'Coincidences', 'Delay')
            #Rt = data.trues_count_rate
            #Rtot = data.prompts_count_rate
            #Rsc = data.scatter_count_rate
            #necr = Rt ** 2 / Rtot
            #sf = Rsc / (Rt + Rsc) * 100
            #Sensitivity = Rt / a * 100
            #print(f'Simulation duration        = {data.duration} sec')
            #print(f'Sensitivity                = {Sensitivity:.2f} %')
            #print(f'ScatterFraction            = {sf:.2f} %   ')
            #print(f'NECR                       = {necr:.2f} cps')
            #activity.append( float(activity_folder[2:-1]) )
            #NECRnoEW.append( necr )
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
            activityEW.append( float(activity_folder[2:-1]) )
            NECREW.append( necr )
        #activity, NECRnoEW = zip(*sorted(zip(activity, NECRnoEW)))
        #activity = np.asarray(activity)
        #NECRnoEW = np.asarray(NECRnoEW)
        activityEW, NECREW = zip(*sorted(zip(activityEW, NECREW)))
        activityEW = np.asarray(activityEW)
        NECREW = np.asarray(NECREW)
        NECREWToF = NECREW*get_NECR_ToFgain(output_folder)
        #print(f'Activity (MBq)  NECR without energy window') 
        #for i in range( n_activities ):
           #print(' %.1f   %.2f ' % (activity[i],NECRnoEW[i]))
        print(f'Activity (MBq)  NECR with energy window (350 - 650 keV)') 
        for i in range( n_activities ):
           print(' %.1f   %.2f ' % (activityEW[i],NECREW[i]))
        #MaxNECRnoEW.append(max(NECRnoEW))
        #MeanNECRnoEW.append(mean(NECRnoEW))
        MaxNECREW.append(max(NECREW))
        MeanNECREW.append(mean(NECREW))
        grNECREW.append(ROOT.TGraph( n_activities, activityEW, NECREW ))
        grNECREW[index].SetLineColor( Color[index] )
        grNECREW[index].SetMarkerColor( Color[index] )
        grNECREW[index].SetTitle(get_legend_text(output_folder))
        grNECR.Add(grNECREW[index])
        grNECREWToF.append(ROOT.TGraph( n_activities, activityEW, NECREWToF ))
        grNECREWToF[index].SetLineColor( Color[index] )
        grNECREWToF[index].SetMarkerColor( Color[index] )
        grNECREWToF[index].SetTitle(get_legend_text(output_folder))
        grNECRToF.Add(grNECREWToF[index])
    print(f' ')
    print("Crystal")
    for i in range( n_folders ):
        print(f'{output_folders[i]},', end=" ")
    print("")
    #print("MaxNECRnoEW")    
    #for i in range( n_folders ):
        #print(f'{MaxNECRnoEW[i]:.2f},', end=" ")
    #print("")
    #print("MeanNECRnoEW")    
    #for i in range( n_folders ):
        #print(f'{MeanNECRnoEW[i]:.2f},', end=" ")
    #print("")
    print("MaxNECREW")    
    for i in range( n_folders ):
        print(f'{MaxNECREW[i]:.2f},', end=" ")
    print("")
    print("MeanNECREW")    
    for i in range( n_folders ):
        print(f'{MeanNECREW[i]:.2f},', end=" ")
    print("")
    c1 = ROOT.TCanvas( 'c1', '', 1700, 1200 )
    c1.SetGridx(1)
    c1.SetGridy(0)
    c1.SetRightMargin(0.02)
    c1.SetTopMargin(0.08)
    c1.SetBottomMargin(0.15)
    c1.SetLeftMargin(0.15)

    pad1 = ROOT.TPad("pad1","pad1" ,0.0 ,0.0 ,0.5 ,0.85)
    #pad1.SetBottomMargin(0)
    pad1.SetLeftMargin(0.08)
    pad1.SetTopMargin(0.04)
    pad1.SetRightMargin(0.02)
    pad1.Draw()
    pad1.cd()
    grNECR.Draw('APC')
    latex1 = ROOT.TLatex()
    latex1.SetNDC()
    latex1.DrawText(0.7 ,0.83 , "non ToF")

    c1.cd()
    pad2 = ROOT.TPad("pad2","pad2" ,0.5 ,0.0 ,1.0 ,0.85)
    #pad2.SetBottomMargin(0)
    pad2.SetLeftMargin(0.09)
    pad2.SetTopMargin(0.04)
    pad2.SetRightMargin(0.02)
    pad2.Draw()
    pad2.cd()
    grNECRToF.Draw('APC')
    latex2 = ROOT.TLatex()
    latex2.SetNDC()
    latex2.DrawText(0.7 ,0.83 , "ToF")

    c1.cd()
    pad3 = ROOT.TPad("pad3","pad3" ,0.01 ,0.86 ,1.0 ,1.0)
    pad3.SetBottomMargin(0.02)
    #pad3.SetTopMargin(0)
    pad3.Draw()
    pad3.cd()
    legendEW = ROOT.TLegend(0., 0., 0.99, 0.99)
    legendEW.SetNColumns(3)
    for index, output_folder in enumerate(output_folders, start=0):
        legendEW.AddEntry(grNECREW[index] ,get_legend_text(output_folder))
    legendEW.Draw('same')

    c1.Update()
    c1.SaveAs(output_folder+"../necr_EWAll.pdf")

    r = True
    return r

def get_NECR_ToFgain(case):
    if case == 'LYSO-Proteus/':
        return 6.35
    elif case == 'BGO/':
        return 2.22
    elif case == 'BGO_BaF2_25/':
        return 5.36
    elif case == 'BGO_EJ232_25/':
        return 6.19
    elif case == 'LYSO-Proteus_BaF2_30/':
        return 11.68
    elif case == 'LYSO-Proteus_EJ232_30/':
        return 11.68
    #elif case == 'BGO_BaF2_23.2/':
    #    return 
    #elif case == 'BGO_EJ232_22.5/':
    #    return 
    #elif case == 'LYSO-Proteus_BaF2_29.7/':
    #    return 
    #elif case == 'LYSO-Proteus_EJ232_31.1/':
    #    return 
    elif case == 'BGO_BaF2_23.7/':
        return 5.39
    elif case == 'BGO_EJ232_22.7/':
        return 6.27
    elif case == 'LYSO-Proteus_BaF2_31.1/':
        return 11.35
    elif case == 'LYSO-Proteus_EJ232_32.3/':
        return 11.02
    else:
        return 1.

def get_legend_text(case):
    if case == 'LYSO-Proteus/':
        return 'LYSO/20 mm'
    elif case == 'BGO/':
        return 'BGO/15 mm'
    elif case == 'BGO_BaF2_25/':
        return 'BGO-BaF2/25 mm'
    elif case == 'BGO_EJ232_25/':
        return 'BGO-EJ232/25 mm'
    elif case == 'LYSO-Proteus_BaF2_30/':
        return 'LYSO-BaF2/30 mm'
    elif case == 'LYSO-Proteus_EJ232_30/':
        return 'LYSO-EJ232/30 mm'
    elif case == 'BGO_BaF2_23.2/':
        return 'BGO-BaF2/23.2 mm'
    elif case == 'BGO_EJ232_22.5/':
        return 'BGO-EJ232/22.5 mm'
    elif case == 'LYSO-Proteus_BaF2_29.7/':
        return 'LYSO-BaF2/29.7 mm'
    elif case == 'LYSO-Proteus_EJ232_31.1/':
        return 'LYSO-EJ232/31.1 mm'
    elif case == 'BGO_BaF2_23.7/':
        return 'BGO-BaF2/23.7 mm'
    elif case == 'BGO_EJ232_22.7/':
        return 'BGO-EJ232/22.7 mm'
    elif case == 'LYSO-Proteus_BaF2_31.1/':
        return 'LYSO-BaF2/31.1 mm'
    elif case == 'LYSO-Proteus_EJ232_32.3/':
        return 'LYSO-EJ232/32.3 mm'
    else:
        return 'no legend'

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
