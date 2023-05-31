# -*- coding: utf-8 -*-
"""
Created on Sun Sep 26 15:29:46 2021

Functions to proccess different deconvoluted file

@author: wingkinlui
"""
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np

#Process topfd msalign file
def topfd(file):
    with open(file) as deconvoluted:
        content = deconvoluted.readlines()
        spectra = {}
        #Iterate lines to get spectrum
        for num, line in enumerate(content):
            if "SCANS" in line:
                scan = int(line.strip("SCANS="))
            if "RETENTION_TIME" in line:
                rt = round(float(line.strip('RETENTION_TIME=').strip('\n'))/60,2)
            if "PRECURSOR_MASS" in line:
                p_mass = float(line.strip("PRECURSOR_MASS=").strip("\n"))
            if "PRECURSOR_INTENSITY" in line:
                p_intensity= float(line.strip("PRECURSOR_INTENSITY=").strip("\n"))
                start_line = num+1
            if "END IONS" in line:
                end_line = num
                
                #Get spectrum
                #Turn spectrum intp np array
                np_spectrum = np.array([i.strip("\n").split("	") for i in content[start_line:end_line]])
                #Convert strings to float
                np_spectrum = np_spectrum.astype(float)
                #Sort by mass
                if len(np_spectrum) >= 2:
                    np_spectrum = np_spectrum[np_spectrum[:,0].argsort()]
                spectrum = {}
                #Combine charge species and put into a dictionary
                #Assign a "memory" i to save the first instance of charge species mass
                memory_i = 0
                for i, fragment in enumerate(np_spectrum):
                    if i == 0:
                        spectrum[fragment[0]] = fragment[1]
                    else:
                        if ((fragment[0] - np_spectrum[i-1][0])*1e6/fragment[0]) <= 5:
                            if memory_i == 0:
                                spectrum[np_spectrum[i-1][0]] += fragment[1]
                                memory_i = i-1
                            else:
                                spectrum[np_spectrum[memory_i][0]] += fragment[1]
                        else:
                            spectrum[fragment[0]] = fragment[1]
                            memory_i = 0
                
                #save scan
                spectra[scan] = {"scan":scan, "rt":rt, "precursor_mass":p_mass, "precursor_intensity":p_intensity,
                                   "spectrum": spectrum}
    deconvoluted.close()
    return spectra

#Process prosight puf file
#Define mass range. Applicable to puf thanks to its xml format

def PC_puf(file, min_mass, max_mass):
    
    tree = ET.parse(file)
    root = tree.getroot()
    spectra = {}
    for exp in root:
        #Get Precursor mass
        p_mass = float(exp[1][2][0][2].text)
        
        
        #skip scan if precursor out of range
        if p_mass < min_mass or p_mass > max_mass:
            continue
        
        #Get precursor intensity
        p_intensity = float(exp[1][2][0][4].text)
        
        #get scan no
        scan = exp[0].text.split("ms2 scans: ")[1].split("]")[0]
        scan = int(scan)
        
        #Get rt
        rt = exp[0].text
        rt = rt.split(") ")[1].split("-")[0]
        rt = float(rt)
        
        
        #Get fragment info
        spectrum = {}
        for fragment in exp[1][3]:
            spectrum[fragment[2].text] = float(fragment[4].text)
        
        #Append scan to spectra dict
        spectra[scan] = {"scan":scan, "rt":rt, "precursor_mass":p_mass, "precursor_intensity":p_intensity,
                                   "spectrum": spectrum}
    return spectra

def PD_puf(file, min_mass, max_mass):
    scan_details = file.replace("_PD.puf", ".raw_Matrix.txt")
    scan_df = pd.read_csv(scan_details, sep = "	")
    
    tree = ET.parse(file)
    root = tree.getroot()
    spectra = {}
    for exp in root:
        #Get Precursor mass
        p_mass = float(exp[1][2][0][3].text)
        
        #skip scan if precursor out of range
        if p_mass < min_mass or p_mass > max_mass:
            continue
        
        #Get precursor intensity
        p_intensity = float(exp[1][2][0][4].text)
        
        #get scan no
        #scan = exp[0].text.split("ms2 scans: ")[1].split("]")[0]
        scan = exp[0].text.split(" ")[1]
        scan = int(scan)
        
        #Get rt
        rt = scan_df[scan_df["MS2ScanNumber"] == scan]["MS2RetTime(min)"].values[0]
        
        #Get fragment info
        spectrum = {}
        for fragment in exp[1][3]:
            spectrum[fragment[2].text] = float(fragment[4].text)
        
        #Append scan to spectra dict
        spectra[scan] = {"scan":scan, "rt":rt, "precursor_mass":p_mass, "precursor_intensity":p_intensity,
                                   "spectrum": spectrum}
    return spectra