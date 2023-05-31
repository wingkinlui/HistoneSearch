# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 11:13:14 2022

This is the script to generate theoretical spectra.
Return a _ms2.msalign file of a theoretical spectra

@author: wingkinlui
"""
import numpy as np


'''
Basic chemical information
'''

#aa masses contains loss of a water molecule. i.e. assume in the middle of a peptide chain -NH...CO-
aa_mass = {
    "A": 71.03711, "R": 156.10111, "N": 114.04293, "D": 115.02694, "C": 103.00919, "E": 129.04259, "Q": 128.05858, "G": 57.02146,
    "H": 137.05891, "I": 113.08406, "L": 113.08406, "K": 128.09496, "M": 131.04049, "F": 147.06841, "P": 97.05276, "S": 87.03203,
    "T": 101.04768, "W": 186.07931, "Y": 163.06333, "V": 99.06841
    }

#C-ion monoisotopic mass: addition of NH3+ from the z-ion and an addition of a H back to the N-term. Minus 1 H for decharge.
c_ion_correction = 17.02654911

#Z-ion monoisotopic mass: addition of an OH back to the C-term, minus a NH since its given to c-ion. 
z_ion_correction = 1.99184063

#Considered amino acids masses
mod_mass = {"unmod": 0, "ac": 42.010565, "me": 14.015650, "me,me": 28.031300, "me,me,me": 42.046950, "ph": 79.966331, "ox2": 31.989829, "ox3": 47.984744}

'''
Parameters
'''
sequence = "ARTKQTARKSTGGKAPRKQLATKAARKSAPSTGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSAAIGALQEASEAYLVGLFEDTNLAAIHAKRVTIAPKDIQLARRIRGERA"

total_mods = {"ac": 1}

proteoforms = [
    [{9: "ac"}, 0.25],
    [{14: "ac"}, 0.25],
    [{23: "ac"}, 0.25],
    [{27: "ac"}, 0.25]
    ]

precursor_intensity = 1000000
fragment_intensity = 10000
deviation = 1e-10  #A deviation range for adding deviations into generated fragment intensity

'''
Create Precursor info
'''
precursor_mass = sum([aa_mass[aa] for aa in sequence]) + 18.0105647 + sum([mod_mass[mod]*total_mods[mod] for mod in total_mods])


'''
Create Spectrum
'''

spectrum = {}

#Iterate proteoform considered
for proteoform in proteoforms:
    modifications = proteoform[0]
    #Iterate sequence length
    ion_dict = {}
    for pos in range(len(sequence)):
        #Get fragment aa composition (including an empty fragment from :0 and len:len)
        #Full length "fragment" is not considered
        c_fragment_aa = sequence[:pos]
        z_fragment_aa = sequence[pos+1:len(sequence)]
        #Get basic ion mass
        c_ion = sum([aa_mass[aa] for aa in c_fragment_aa]) + c_ion_correction
        z_ion = sum([aa_mass[aa] for aa in z_fragment_aa]) + z_ion_correction
        for mod_pos in modifications:
            #Add mod mass if mod is within c-ion
            if mod_pos <= pos:
                c_ion += mod_mass[modifications[mod_pos]]
            #Add mod mass if mod is within z-ion
            if mod_pos > pos+1:
                z_ion += mod_mass[modifications[mod_pos]]
        #Other than the empty aa mass, add masses and corresponding intenities into a temp ion dictionary
        if c_ion != c_ion_correction:
            ion_dict[c_ion] = fragment_intensity*proteoform[1] + np.random.uniform(0, fragment_intensity*deviation)
        if z_ion != z_ion_correction:
            ion_dict[z_ion] = fragment_intensity*proteoform[1] + np.random.uniform(0, fragment_intensity*deviation)
    #Add ions into spectrum
    for ion in ion_dict:
        if ion not in spectrum:
            spectrum[ion] = ion_dict[ion]
        else:
            spectrum[ion] += ion_dict[ion]
        
'''
Build msalign file
'''        

name = str(round(precursor_mass))
for mod in total_mods:
    name += "_" + str(total_mods[mod])+mod

header = open("msalign_header.txt","r")
msalign = open(name + "_ms2.msalign", "w")
for line in header.readlines():
    msalign.write(line)

msalign.write("PRECURSOR_MZ="+str(round(precursor_mass/20, 2))+"\n")
msalign.write("PRECURSOR_CHARGE=20"+"\n")
msalign.write("PRECURSOR_MASS="+str(precursor_mass)+"\n")
msalign.write("PRECURSOR_INTENSITY="+str(precursor_intensity)+"\n")

for frag in spectrum:
    msalign.write(str(frag)+"	"+str(spectrum[frag])+"	"+str(1)+"\n")

msalign.write("END IONS"+"\n\n")

msalign.close()



    
