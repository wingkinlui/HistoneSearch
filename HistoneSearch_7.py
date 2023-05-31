# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 12:00:39 2021

Main function to load deconvoluted spectra, FASTA, and other modules.
Automatically read .puf or _ms2.msalign file, and FASTA files in the directory.
Searched spectra and their PrSMs would be stored as .pkl file. If they exist in the directory, the .pkl file would be read and directly go to quantification

variable PTMs and fixed modification to be specified in the Search_Function module.

Return searched spectra (.pkl), identified and quantified proteoforms for each MS/MS (.csv), and error report file (.csv)

@author: wingkinlui
"""


import os
import pandas as pd
import csv
import pickle

#from Quantifier_function_2 import Quantifer
from Search_Functions_6_1 import search, mod_mass, fixed_mods, aa_mass, t_tol
from Deconvolution_File_Processor_2 import topfd, PC_puf, PD_puf
from FIRR_Quantification_1 import Quantification
'''
Filter section
'''
#Precursor mass milter
min_mass = 11000
max_mass = 12000

min_rt = 0
max_rt = 180


'''
Run section
'''

files = []
#Read files in directory to get deconvoluted files, searched file and fasta
for file in os.listdir("."):
    if file.endswith(".puf") or file.endswith("_ms2.msalign") or file.endswith(".pkl"):
        files.append(file)        
    if file.endswith(".fasta"):
        with open(file) as fasta:
            fasta_content = fasta.readlines()
            for num, line in enumerate(fasta_content):
                if ">" in line:
                    sequence = fasta_content[num+1].strip("M")

#After getting fasta, get basic mass and total fixed mod mass
#The aa mass list does not contain the water molecule mass. Its mass will be added here.
basic_mass = 1.007825 + 1.007825 + 15.994915
for resid in sequence:
    basic_mass += aa_mass[resid]

fixed_mod_mass = 0
for resid in sequence:
    if resid in fixed_mods:
        fixed_mod_mass += mod_mass[fixed_mods[resid]]

'''
Deconvolution file processing
'''

#Iterate each files to get experiments
experiments = []
for file in files:
    #Get names of the files
    if file.endswith("_ms2.msalign"):
        name = file.replace("_ms2.msalign", "")
    if file.endswith(".puf"):
        name = file.replace(".puf", "")
    if file.endswith(".pkl"):
        name = file.replace(".pkl", "")
    experiments.append(name)
        
experiments = set(experiments)

#Iterate experiment
for name in experiments:
    #Skip search if searched file already present
    if name + ".pkl" not in files:
        #Read deconvoluted file and get spectra
        if name + "_ms2.msalign" in files:
            spectra = topfd(name + "_ms2.msalign")
        if name + ".puf" in files:
            spectra = PC_puf(name + ".puf", min_mass, max_mass)
            
        search_results = {}
        
        #Create a progress variable to show progress    
        Progress = 0
        total_scans = len(spectra)
        
        '''
        Proteoform Searching
        '''
        
        for scan in spectra:
            matched_candidates = []
            candidate = {
                    "cur_pos" : 0,
                    "cur_th_mass" : 0,
                    "unresolved" : {},
                    "modifications" : {},
                    "matches" : {},
                    "Continue" : True,
                    "ac_mod" : 0,
                    "ac_mod_mass": 0
                    }
            #mass filter and perform search; save mod states for post-search hit processing
            if spectra[scan]["precursor_mass"] > min_mass and spectra[scan]["precursor_mass"] < max_mass:
                if spectra[scan]["rt"] >= min_rt and spectra[scan]["rt"] <= max_rt:
                    max_mod_mass = spectra[scan]["precursor_mass"] - basic_mass - fixed_mod_mass + t_tol
                    #if max_mod_mass <=0:
                        #continue
                    spectrum = spectra[scan]["spectrum"]
                    search(candidate, spectrum, sequence, max_mod_mass, basic_mass, fixed_mod_mass, spectra[scan]["precursor_mass"], matched_candidates)
                    #Pass scan if there is no matched candidates
                    if len(matched_candidates) == 0:
                        continue
                    else:
                        search_results[scan] = spectra[scan]
                        search_results[scan]["candidates"] = matched_candidates
                
            Progress += 1
            print("Searching. Current Scan " + str(scan) + '; '+ str(Progress) + " out of " + str(total_scans) + " scans completed.")
    
        
        '''
        Save search results
        '''
        dict_file = open(name+".pkl", "wb")
        pickle.dump(search_results, dict_file)
        dict_file.close()
    
    
    '''
    Post Search Processing and Quantification
    '''
    
    if name + ".pkl" in files: 
        dict_file = open(name + ".pkl", "rb")
        search_results = pickle.load(dict_file)
        dict_file.close()
    
    #Create a report dictionary for quantification report
    report = {"Quantified": [], "NoInformativeGrp" : [], "UndefinedSol": [], "Converged": [], "All_Filtered": []}
   
    print("Post-search processing and quantification started")
    
    post_process_results, report = Quantification(search_results, report)
    
    '''
    #Skip co-isolated qaunt
    post_process_results = search_results
    for scan in search_results:
        #search_results[scan]["candidates"] = [search_results[scan]["candidates"][0]]
        for candidate_i in range(len(search_results[scan]["candidates"])):
            search_results[scan]["candidates"][candidate_i]["Abundance"] = search_results[scan]["precursor_intensity"]
    '''
    
    #Drop scans that has no candidate after post search processing
    
    print("Post-search processing and quantification finished")
    
    #Convert quantified results to dataframe
    print("Writing results")
    quanti_result = []
    ID = 1
    for scan in post_process_results:
        for candidate in post_process_results[scan]["candidates"]:
            quanti_result.append([ID, candidate["name"], scan, candidate["final_th_mass"], candidate["Abundance"], post_process_results[scan]["rt"], len(candidate["matches"]), list(candidate["matches"].keys())[-1] - list(candidate["matches"].keys())[0]])
            ID += 1
    Summarydf = pd.DataFrame(quanti_result, columns = ["ID", "Proteoform", "Scan", "Mass", "Abundance", "RT", "#MatchedIons", "MatchedIonRange"])
    
    
    '''
    Writing up results
    '''
    
    #Convert report to datafram
    for key in report:
        report[key] = list(set(report[key]))
    with open(name + "_Report.csv", "w") as report_csv:
        writer = csv.writer(report_csv)
        for key, value in report.items():
            row = [key] + [i for i in value]
            writer.writerow(row)
    report_csv.close()
    
    
    #Export df to csv
    Summarydf.to_csv(name + "_Quantified.csv", index = False)
    
    print("Finished\n")
         