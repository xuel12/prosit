#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 20:51:37 2020

@author: xuel12
"""

import os
os.chdir('/Users/xuel12/Documents/MSdatascience/CS7180AI/project/prosit/')
import re


def from_msp_prosit(ofile):
    
    # Get a list of the keys for the datasets
    f = open('examples/peptidelist.msp', 'r')
    tmp_dict = {}
    final_dict = {}
    # use readline() to read the first line 
    line = f.readline()
    # use the read line to read further.
    # If the file is not empty keep reading one line at a time, till the file is empty
    while line:
        if re.search("^Name",line) is not None:
            tmp_dict['Name'] = line.rstrip('\n').split(' ')[-1]
            line = f.readline()
        elif re.search("^MW",line) is not None:
            tmp_dict['MW'] = round(float(line.rstrip('\n').split(' ')[-1]),4)
            line = f.readline()
        elif re.search("^Comment",line) is not None:
            tmp_dict['Comment'] = re.sub('Comment: ','',line.rstrip('\n'))
            tmpline = tmp_dict['Comment'].split(' ')
            parent, ce, mod, modseq = [x.split('=')[1] for x in tmpline]
            ce = round(float(ce),1)
            line = f.readline()
        elif re.search("^Num peaks",line) is not None:
            tmp_dict['Num_peaks'] = line.rstrip('\n').split(' ')[-1]
            line = f.readline()
            mz = []
            intensity = []
            anno = []
            while (line and line.strip() and re.search("^Name",line) is None):
                tmpline = line.rstrip('\n').split('\t')
                mz.append(tmpline[0])
                intensity.append(tmpline[1])
                anno.append(tmpline[2])
                line = f.readline()
            # one record per peptide_ce_modstring
            peptide_ce_modstring = '_' .join([tmp_dict['Name'],str(tmp_dict['MW']), str(ce), mod])
#            final_dict['peptide_ce_modstring'] = '_' .join([tmp_dict['Name'],tmp_dict['Comment'][1],tmp_dict['Comment'][3]])
            final_dict[peptide_ce_modstring] = {}
            final_dict[peptide_ce_modstring]['mz'] = mz
            final_dict[peptide_ce_modstring]['intensity'] = intensity
            final_dict[peptide_ce_modstring]['anno'] = anno
    f.close()
    return final_dict
    
spectrum_prosit = from_msp_prosit('examples/peptidelist.msp')

def from_msp_propel(ofile):
    
    # Get a list of the keys for the datasets
    with open(ofile, 'r') as f:
        tmp_dict = dict()
        final_dict = {}
        # use readline() to read the first line 
        line = f.readline()
        # use the read line to read further.
        # If the file is not empty keep reading one line at a time, till the file is empty
        counter = 0
        max_count = 100
        ce = 0
        while (line is not None):
            if not line.strip():
                line = f.readline()
            elif counter >= max_count:
                break
            else:
                if re.search("^Name",line) is not None:
                    tmpline = line.rstrip('\n').split(' ')[-1]
                    tmp_dict['Name'] = tmpline.split('_')[0]
                    # print(tmp_dict['Name'])
                    mod, ce = [tmpline.split('_')[1], tmpline.split('_')[2]]
                    ce = round(float(ce.strip('%')),1)
                    line = f.readline()
                elif re.search("^MW",line) is not None:
                    tmp_dict['MW'] = round(float(line.rstrip('\n').split(' ')[-1]),4)
                    line = f.readline()
                elif re.search("^Comment",line) is not None:
                    tmpline = re.sub('Comment: ','',line.rstrip('\n'))
                    # tmpline = line.rstrip('\n').split(' ')[1:]
                    # tmpline = [x.split('=')[1] for x in tmpline]
                    tmp_dict['Comment'] = tmpline
                    line = f.readline()
                elif re.search("^Num peaks",line) is not None:
                    tmp_dict['Num_peaks'] = line.rstrip('\n').split(' ')[-1]
                    line = f.readline()
                    mz = []
                    intensity = []
                    anno = []
                    while (line is not None):
                        if not line.strip():
                            # one record per peptide_ce_modstring
                            peptide_ce_modstring = '_'.join([tmp_dict['Name'],str(tmp_dict['MW']), str(ce), mod])
                            final_dict[peptide_ce_modstring] = {}
                            final_dict[peptide_ce_modstring]['mz'] = mz
                            final_dict[peptide_ce_modstring]['intensity'] = intensity
                            final_dict[peptide_ce_modstring]['anno'] = anno
                            # tmp_dict = dict()
                            counter += 1
                            break
                        else:
                            tmpline = line.rstrip('\n').split('\t')
                            mz.append(tmpline[0])
                            intensity.append(tmpline[1])
                            anno.append(tmpline[2])
                            line = f.readline()
    return final_dict

spectrum_propel = from_msp_propel('examples/human_synthetic_hcd_selected.msp')


def to_mgf(spectrum_dict, ofile):   
    # take spedtrum list, write to mgf
    with open(ofile, 'w') as f:
        rtinseconds = 1
        for name in spectrum_dict:
            f.write('BEGIN IONS\n')
            tmpname = name.split('_')
            seq = tmpname[0].split('/')[0]
            charge = tmpname[0].split('/')[1]
            pepmass = tmpname[1]
            f.write('SEQ={}\n'.format(seq))
            f.write('PEPMASS={}\n'.format(pepmass))
            f.write('CHARGE={}+\n'.format(charge))
            f.write('TITLE={}\n'.format(name))
            f.write('RTINSECONDS={}\n'.format(rtinseconds))
            for i in range(len(spectrum_dict[name]['mz'])):
                tmpion = spectrum_dict[name]['mz'][i] + ' ' + spectrum_dict[name]['intensity'][i] + '\n'
                f.write(tmpion)
            f.write('END IONS\n')
            rtinseconds += 1
    f.close()
    
to_mgf(spectrum_dict, 'examples/peptidelist.mgf')


