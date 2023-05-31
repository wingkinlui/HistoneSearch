# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 11:13:14 2022

This is the script to generate theoretical spectra.
Return a .puf file.

@author: wingk
"""
import numpy as np
import xml.etree.ElementTree as ET
from xml.dom import minidom

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
sequence = "GARTKQTARKSTGGKAPRKQLATKAARKSAPSTGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSAAIGALQEASEAYLVGLFEDTNLAAIHAKRVTIAPKDIQLARRIRGERA"

total_mods = {"ac": 2}

proteoforms = [
    [{5:"ac", 10: "ac"}, 0.5],
    [{24:"ac", 37: "ac"}, 0.5]
    ]

precursor_intensity = 1000000
fragment_intensity = 10000
deviation = 0.01  #A deviation range for adding deviations into generated fragment intensity

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
Build XML (puf) file
'''        
new_root = ET.Element("data_set")
new_root.set("owner", "")
new_root.set("type", "")
new_root.set("version", "1.1")

MS_Exp = ET.SubElement(new_root, "ms-ms_experiment", attrib={"id":"0", "source":""})

comment = ET.SubElement(MS_Exp, "comment")
comment.text = "ETD fragmentation for precursor at m/z " + str(round(precursor_mass/20, 2)) + " from retention time (min) 0-0.1 [ms1 scans: 998,1000; ms2 scans: 999] with FourierTransform detection."

instrument_data = ET.SubElement(MS_Exp, "instrument_data")

frag_method = ET.SubElement(instrument_data, "fragmentation_method")
frag_method.text = "ETD"

ion_type = ET.SubElement(instrument_data, "ion_type")
ion_type.text = "CZ"

intact_list = ET.SubElement(instrument_data, "intact_list")

intact = ET.SubElement(intact_list, "intact", attrib={"id":"1"})

mz_mono = ET.SubElement(intact, "mz_monoisotopic")
mz_mono.text = str(precursor_mass/20)

mz_aver = ET.SubElement(intact, "mz_average")
mz_aver.text = "0"

mass_mono = ET.SubElement(intact, "mass_monoisotopic")
mass_mono.text = str(precursor_mass)

mass_aver = ET.SubElement(intact, "mass_average")
mass_aver.text = "0"

mass_int = ET.SubElement(intact, "intensity")
mass_int.text = str(precursor_intensity)

frag_list = ET.SubElement(instrument_data, "fragment_list")

id_no = 1

for fragment in spectrum:
    frag = ET.SubElement(frag_list, "fragment", attrib={"id":str(id_no)})
    
    frag_mz_mono = ET.SubElement(frag, "mz_monoisotopic")
    frag_mz_mono.text = "0"
    
    frag_mz_aver = ET.SubElement(frag, "mz_average")
    frag_mz_aver.text = "0"
    
    frag_mass_mono = ET.SubElement(frag, "mass_monoisotopic")
    frag_mass_mono.text = str(fragment)
    
    frag_mass_aver = ET.SubElement(frag, "mass_average")
    frag_mass_aver.text = "0"
    
    frag_int = ET.SubElement(frag, "intensity")
    frag_int.text = str(spectrum[fragment])
    
    id_no += 1

analysis = ET.SubElement(MS_Exp, "analysis", attrib = {"id": "0", "type": "null"})
analysis_comment = ET.SubElement(analysis, "comment")
analysis_comment.text = "Null Search"

#Pretty print
xmlstr = minidom.parseString(ET.tostring(new_root)).toprettyxml(indent="   ")

#Write file
name = str(round(precursor_mass))
for mod in total_mods:
    name += "_" + str(total_mods[mod])+mod
    
with open(name + ".puf", "w") as f:
    f.write(xmlstr)

    
