#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 22:50:49 2020

@author: xuel12
"""
import os
os.chdir('/Users/xuel12/Documents/MSdatascience/CS7180AI/project/prosit/local_training')
import matplotlib
import matplotlib.pyplot as plt
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
from pyteomics import mgf
import numpy as np
import constants
import msp_parser

##### single plot
def singleplot(feature, file):
    # Read the spectrum from an MGF file using Pyteomics.
    spectrum_dict = mgf.get_spectrum(file, feature)
    # modifications = {8: 15.994915}
    modifications = {}
    
    identifier = spectrum_dict['params']['title']
    precursor_mz = spectrum_dict['params']['pepmass'][0]
    precursor_charge = spectrum_dict['params']['charge'][0]
    mz = spectrum_dict['m/z array']
    intensity = spectrum_dict['intensity array']
    retention_time = float(spectrum_dict['params']['rtinseconds'])
    peptide = spectrum_dict['params']['seq']
    
    # # if sorted mz is desired
    # pair = np.vstack([spectrum_dict['m/z array'], spectrum_dict['intensity array']]).T
    # sorted_pair = pair[np.argsort(pair[:, 0])]
    # mz, intensity = [sorted_pair[:,0], sorted_pair[:,1]]
    
    # Create the MS/MS spectrum.
    spectrum = sus.MsmsSpectrum(identifier, precursor_mz, precursor_charge, mz, intensity, \
                                retention_time=retention_time, peptide=peptide, \
                                modifications=modifications)
    
    # Process the MS/MS spectrum.
    fragment_tol_mass = 0.5
    fragment_tol_mode = 'Da'
    min_mz = 100
    min_intensity = 0.05
    spectrum = (spectrum.set_mz_range(min_mz=min_mz, max_mz=1400)
                .remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
                .filter_intensity(min_intensity=min_intensity, max_num_peaks=50)
                # .scale_intensity('root')
                .annotate_peptide_fragments(fragment_tol_mass, fragment_tol_mode, ion_types='aby'))
    # # label the mz for all peaks
    # annotate_fragment_mz = sorted_pair[(sorted_pair[:,0]>min_mz) & (sorted_pair[:,1]>min_intensity), 0]
    # for fragment_mz in annotate_fragment_mz:
    #     spectrum.annotate_mz_fragment(fragment_mz, 1, fragment_tol_mass, fragment_tol_mode)
        
    # Plot the MS/MS spectrum.
    fig, ax = plt.subplots(figsize=(12, 6))
    plt.title(identifier)
    sup.spectrum(spectrum, ax=ax)
    plt.show()
    plt.close()


##### mirror plot
def mirroplot_twopeptides(peplist, file):
    fragment_tol_mass = 0.5
    fragment_tol_mode = 'Da'
    spectra = []
    for spectrum_dict in mgf.read(file):
        if peplist[0] in spectrum_dict['params']['title'] or peplist[1] in spectrum_dict['params']['title']:
            identifier = spectrum_dict['params']['title']
            precursor_mz = spectrum_dict['params']['pepmass'][0]
            precursor_charge = spectrum_dict['params']['charge'][0]
            mz = spectrum_dict['m/z array']
            intensity = spectrum_dict['intensity array']
            retention_time = float(spectrum_dict['params']['rtinseconds'])
            peptide = spectrum_dict['params']['seq']
            # modifications = {6: 15.994915}
            modifications = {}
    
            # Create the MS/MS spectrum.
            spectra.append(sus.MsmsSpectrum(identifier, precursor_mz,
                                            precursor_charge, mz, intensity,
                                            retention_time=retention_time,
                                            peptide=peptide,
                                            modifications=modifications)
                           .filter_intensity(0.01, 50)
                           # .scale_intensity('root')
                           .annotate_peptide_fragments(fragment_tol_mass, fragment_tol_mode, ion_types='aby'))
                
    fig, ax = plt.subplots(figsize=(12, 6))
    spectrum_top, spectrum_bottom = spectra
    sup.mirror(spectrum_top, spectrum_bottom, ax=ax)
    plt.show()
    plt.close()


##### mirror plot for two dataset
def mirroplot_twosets(peplist, pred_file, ref_file, plot_dir):
    import re
    
    min_intensity = 0.01
    fragment_tol_mass = 0.5
    fragment_tol_mode = 'Da'
    title = peplist[4]
    for title in peplist:
        spectra = []
        pred_dict = mgf.get_spectrum(example_dir+'peptidelist_pred.mgf', title)
        ref_dict = mgf.get_spectrum(data_dir+'human_synthetic_hcd_selected.mgf', title)
        if (ref_dict is None):
            break
        pair = [pred_dict, ref_dict]
        for spectrum_dict in pair:
            identifier = spectrum_dict['params']['title']
            precursor_mz = spectrum_dict['params']['pepmass'][0]
            precursor_charge = spectrum_dict['params']['charge'][0]
            mz = spectrum_dict['m/z array']
            intensity = spectrum_dict['intensity array']
            retention_time = float(spectrum_dict['params']['rtinseconds'])
            peptide = spectrum_dict['params']['seq']
            # modifications = {6: 15.994915}
            modifications = {}
            
            # Create the MS/MS spectrum.
            spectra.append(sus.MsmsSpectrum(identifier, precursor_mz,
                                            precursor_charge, mz, intensity,
                                            retention_time=retention_time,
                                            peptide=peptide,
                                            modifications=modifications)
                           .filter_intensity(min_intensity, 50)
                           # .scale_intensity('root')
                           .annotate_peptide_fragments(fragment_tol_mass, fragment_tol_mode, ion_types='aby'))
            
        fig, ax = plt.subplots(figsize=(12, 6))
        plt.title(identifier)
        spectrum_top, spectrum_bottom = spectra
        sup.mirror(spectrum_top, spectrum_bottom, ax=ax)
        fig.savefig(plot_dir+'/{}.png'.format(re.sub('/','_',identifier)))
        plt.close(fig)
        # plt.show()
        # plt.close()


def peplist_from_csv(csvfile):
    peptidelist = []
    with open (csvfile, 'r') as f:
        f.readline()
        for line in f:
            seq, ce, charge = line.rstrip('\n').split(',')
            peptide = seq + '/' + charge + '_' + str(round(float(ce),1)) + '_' + '0'
            peptidelist.append(peptide)
    return (peptidelist)
    
    
    
if __name__ == "__main__":
    os.chdir(constants.BASE_PATH + 'project/prosit/local_training')
    data_dir = constants.DATA_DIR
    example_dir = constants.EXAMPLE_DIR
    plot_dir = constants.PLOT_DIR
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    
    # get list of peptides for plotting
    peplist = peplist_from_csv(example_dir + '/peptidelist.csv')

    # store msp files to dictionary from prosit prediction
    spectrum_prosit = msp_parser.from_msp_prosit(example_dir+'peptidelist_pred.msp')
    msp_parser.dict2mgf(spectrum_prosit, example_dir+'peptidelist_pred.mgf')
    
    # single spectra
    singleplot(peplist[0], example_dir+'peptidelist_pred.mgf')
    # compare two different peptides
    mirroplot_twopeptides(peplist[:2], example_dir+'peptidelist_pred.mgf')
    # compare same peptide from two methods
    mirroplot_twosets(peplist, example_dir+'peptidelist_pred.mgf', data_dir+'human_synthetic_hcd_selected.mgf', plot_dir)
        

        