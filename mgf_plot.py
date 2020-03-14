#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 22:50:49 2020

@author: xuel12
"""
import os
os.chdir('/Users/xuel12/Documents/MSdatascience/CS7180AI/project/prosit/')
import matplotlib.pyplot as plt
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
from pyteomics import mgf
import numpy as np

##### single plot

# Read the spectrum from an MGF file using Pyteomics.
spectrum_dict = mgf.get_spectrum('examples/peptidelist.mgf', \
                                 'MLAPPPIMK/2_507.27974918115007_30.0_MLAPPPIMK//Oxidation@M8/2')
modifications = {8: 15.994915}

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
sup.spectrum(spectrum, ax=ax)
plt.show()
plt.close()


##### mirror plot
spectra = []
for spectrum_dict in mgf.read('examples/peptidelist.mgf'):
    if 'MLAPPPIMK' in spectrum_dict['params']['seq'] or 'MRALLLIPPPPMR' in spectrum_dict['params']['seq']:
        identifier = spectrum_dict['params']['title']
        precursor_mz = spectrum_dict['params']['pepmass'][0]
        precursor_charge = spectrum_dict['params']['charge'][0]
        mz = spectrum_dict['m/z array']
        intensity = spectrum_dict['intensity array']
        retention_time = float(spectrum_dict['params']['rtinseconds'])
        peptide = spectrum_dict['params']['seq']
        modifications = {6: 15.994915}

        # Create the MS/MS spectrum.
        spectra.append(sus.MsmsSpectrum(identifier, precursor_mz,
                                        precursor_charge, mz, intensity,
                                        retention_time=retention_time,
                                        peptide=peptide,
                                        modifications=modifications)
                       .filter_intensity(0.01, 50)
                       # .scale_intensity('root')
                       .annotate_peptide_fragments(0.5, 'Da', ion_types='aby'))
        
        
fig, ax = plt.subplots(figsize=(12, 6))
spectrum_top, spectrum_bottom = spectra
sup.mirror(spectrum_top, spectrum_bottom, ax=ax)
plt.show()
plt.close()