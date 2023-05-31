# -*- coding: utf-8 -*-
"""
Created on Sat Sep 25 16:29:17 2021

Functions to perform proteoform search

Variable PTMs and fixed modification to be specified here.
Also specify the mass tolerance, allowed number of modifications, coverage tolerance, etc.

@author: wingkinlui
"""
import copy

'''
Search parameters
'''

f_tol = 5 #In ppm
t_tol = 2 #In Da
max_no_mod = 7 #Max. no of mod to occur on a proteoform
min_frag = 15 #Hit with fragment below this number will be filtered out
mod_coverage_tol = 10 #Hit with mod fragment gap above this range will be discarded
strict_coverage_range = 0 #Under this coverage range, even unmod modifiable residue would also require site determination ion



#Dictionary of amino acid masses, N1H1 at N and O1 at C
aa_mass = {
    "A": 71.03711, "R": 156.10111, "N": 114.04293, "D": 115.02694, "C": 103.00919, "E": 129.04259, "Q": 128.05858, "G": 57.02146,
    "H": 137.05891, "I": 113.08406, "L": 113.08406, "K": 128.09496, "M": 131.04049, "F": 147.06841, "P": 97.05276, "S": 87.03203,
    "T": 101.04768, "W": 186.07931, "Y": 163.06333, "V": 99.06841
    }

#Considered amino acids masses
mod_mass = {"unmod": 0, "ac": 42.010565, "me": 14.015650, "me,me": 28.031300, "me,me,me": 42.046950, "ph": 79.966331, "ox2": 31.989829, "ox3": 47.984744}
#Considered dynamic mods
#REMEMBER TO PUT "UNMOD"!!!!!
#residue_mods = {"K": ["ac", "unmod", "me", "me,me", "me,me,me"], "R": ["unmod", "me", "me,me"], "S": ["unmod","ph"], "T": ["unmod","ph"]}
residue_mods = {"K": ["ac", "unmod", "me", "me,me", "me,me,me"]}

#N-term mods
n_mods = ["ac", "unmod"]

#Fixed mods
fixed_mods = {"M": "ox2", "C": "ox3"}

'''
Search Functions
'''
#Get minimium possible mod mass for capping in later step
min_mod_mass = []
for i in mod_mass:
    if i not in fixed_mods.values() and i != "unmod":
        min_mod_mass.append(mod_mass[i])
min_mod_mass = min(min_mod_mass)


#This is a function to try match current aa mass + next aa mass + putative mod mass
#If there is no match, NoneType will be return. Otherwise, matched fragment intensity will be returned.
#Candidate dict will not be changed here
def match_ion(species, candidate, aa, mod, spectrum):
    #get current aa
    #Update theoretical mass
    attempt = candidate["cur_th_mass"] + aa_mass[aa] + mod_mass[mod]
    #print(attempt)
    match = None
    #Perform mass match
    for mass in spectrum.keys():
        flt_mass = float(mass)
        if abs((attempt+species-flt_mass)*1000000/(attempt+species)) < f_tol:
            match = (mass, spectrum[mass])
            break;
    return match

#This is a function to take match value from match_ion and change candidate dict base on match and aa
def decision(match, candidate, aa, mod, max_mod_mass):
    
    #Create new candidate dict to avoid overwriting
    new_candidate = copy.deepcopy(candidate)
    
    #If current aa has fixed mod:
    if aa in fixed_mods:
        #If matched, add matched fragment to matches dict
        if match != None:
            new_candidate["matches"][candidate["cur_pos"]+1] = [aa, match]
            new_candidate["cur_pos"] += 1
            new_candidate["cur_th_mass"] += aa_mass[aa] + mod_mass[mod]
            return new_candidate
        #If not matched, update position values
        else:
            new_candidate["cur_pos"] += 1
            new_candidate["cur_th_mass"] += aa_mass[aa] + mod_mass[mod]
            return new_candidate
        
    #If current aa could be modified
    elif aa in residue_mods or (candidate["cur_pos"] == 0 and len(n_mods) > 0):
        #print(candidate)
        #If pos within strict coverage range and previously unresolved, cut no matter what
        if candidate["cur_pos"] <= strict_coverage_range and candidate["unresolved"] != {}:
            new_candidate["Continue"] = False
            return new_candidate
        
        #If previously unresolved but position is out of strict coverage range,
        elif candidate["unresolved"] != {}:
            
            #if previously looking for mod
            if candidate["unresolved"]["mod"] != "unmod":
                #If previously looking for mod and currently attempting a mod, cut
                if mod != "unmod":
                    new_candidate["Continue"] = False
                    return new_candidate
                
                #If previously looking for mod and the mod can occur in current aa, cut
                elif candidate["unresolved"]["mod"] in residue_mods[aa]:
                    new_candidate["Continue"] = False
                    return new_candidate
                
                #if previously looking for mod and currently not search for mod, attempt match:
                elif mod == "unmod":
                    #If matched
                    if match != None:
                        #append modified matched info
                        new_candidate["matches"][candidate["cur_pos"]+1] = [aa, match]
                        new_candidate["modifications"][candidate["unresolved"]["pos"]+1] = [candidate["unresolved"]["aa"], candidate["unresolved"]["mod"]]
                        #Return unresolved to {}
                        new_candidate["unresolved"] = {}
                        new_candidate["cur_pos"] += 1
                        new_candidate["cur_th_mass"] += aa_mass[aa]
                        return new_candidate
                    #If not matched, unresolved will remain to be the first instance of mod attempt
                    else:
                        new_candidate["unresolved"]["span"] += 1
                        new_candidate["cur_pos"] += 1
                        new_candidate["cur_th_mass"] += aa_mass[aa]
                        return new_candidate
                
            
            #If previously not looking for mod and currently search for mod,
            elif candidate["unresolved"]["mod"] == "unmod" and mod != "unmod":     
                
                #Get possible mods from unresolved
                possible_mods = candidate["unresolved"]["possible_mods"]
                        
                #currently matching a mod that could also occur in unresolved aa, cut
                if mod in possible_mods:
                    new_candidate["Continue"] = False
                    #print(new_candidate)
                    return new_candidate
                
                #if current mod not possible on previous unresolved aa, attempt match
                elif mod not in possible_mods:
                    #If matched
                    if match != None:
                        #append modified matched info
                        new_candidate["matches"][candidate["cur_pos"]+1] = [aa, match]
                        new_candidate["modifications"][candidate["cur_pos"]+1] = [aa, mod]
                        #Return unresolved to {}
                        new_candidate["unresolved"] = {}
                        new_candidate["ac_mod"] += 1
                        new_candidate["ac_mod_mass"] += mod_mass[mod]
                        new_candidate["cur_pos"] += 1
                        new_candidate["cur_th_mass"] += aa_mass[aa] + mod_mass[mod]
                        return new_candidate
                        
                    #If not matched, update unresolved
                    else:
                        new_candidate["unresolved"]["aa"] = aa
                        new_candidate["unresolved"]["pos"] = candidate["cur_pos"]
                        new_candidate["unresolved"]["mod"] = mod
                        new_candidate["unresolved"]["span"] += 1
                        new_candidate["ac_mod"] += 1
                        new_candidate["ac_mod_mass"] += mod_mass[mod]
                        new_candidate["cur_pos"] += 1
                        new_candidate["cur_th_mass"] += aa_mass[aa] + mod_mass[mod]
                        return new_candidate
                
            #else, (if unresolved but unmod and looking for unmod):
            else:
                #If matched
                if match != None:
                    new_candidate["matches"][candidate["cur_pos"]+1] = [aa, match]
                    new_candidate["cur_pos"] += 1
                    new_candidate["cur_th_mass"] += aa_mass[aa]
                    new_candidate["unresolved"] = {}
                    return new_candidate
                
                #If not matched, unresolved will remain to be the first instance, possible mods will be updated
                else:
                    new_candidate["unresolved"]["span"] += 1
                    new_candidate["unresolved"]["pos"] = candidate["cur_pos"]
                    new_candidate["unresolved"]["mod"] = mod
                    new_candidate["unresolved"]["possible_mods"] += [i for i in residue_mods[aa] if i != "unmod"]
                    new_candidate["cur_pos"] += 1
                    new_candidate["cur_th_mass"] += aa_mass[aa]
                    return new_candidate
            
        #If previously resolved:
        else:
            #If matched
            if match != None:
                #If matching to modification, append modified matched info
                if mod != "unmod":
                    new_candidate["matches"][candidate["cur_pos"]+1] = [aa, match]
                    new_candidate["modifications"][candidate["cur_pos"]+1] = [aa, mod]
                    #Return unresolved to ()
                    new_candidate["unresolved"] = {}
                    new_candidate["ac_mod"] += 1
                    new_candidate["ac_mod_mass"] += mod_mass[mod]
                    new_candidate["cur_pos"] += 1
                    new_candidate["cur_th_mass"] += aa_mass[aa] + mod_mass[mod]
                    return new_candidate
                #If matching to unmod, append unmodified matched info
                else:
                    new_candidate["matches"][candidate["cur_pos"]+1] = [aa, match]
                    new_candidate["cur_pos"] += 1
                    new_candidate["cur_th_mass"] += aa_mass[aa]
                    new_candidate["unresolved"] = {}
                    return new_candidate
                
            #If nothing is matched add position to unresolved and update positions
            else:
                #if attempting mod
                if mod != "unmod":
                    new_candidate["unresolved"] = {"aa": aa, "pos": candidate["cur_pos"], "mod": mod, "span": 1}
                    new_candidate["ac_mod"] += 1
                    new_candidate["ac_mod_mass"] += mod_mass[mod]
                    new_candidate["cur_pos"] += 1
                    new_candidate["cur_th_mass"] += aa_mass[aa] + mod_mass[mod]
                    return new_candidate
                
                #if attempting unmod
                else:
                    #Special handling if currently search for n term mod
                    if candidate["cur_pos"] == 0:
                        new_candidate["unresolved"] = {"aa": aa, "pos": candidate["cur_pos"], "mod": mod, "possible_mods": [i for i in n_mods if i != "unmod"], "span": 1}
                    else:
                        new_candidate["unresolved"] = {"aa": aa, "pos": candidate["cur_pos"], "mod": mod, "possible_mods": [i for i in residue_mods[aa] if i != "unmod"], "span": 1}
                    new_candidate["cur_pos"] += 1
                    new_candidate["cur_th_mass"] += aa_mass[aa]
                    return new_candidate
        
        
    #if aa could not be modified but there is previously unresolved modifible residue
    elif aa not in residue_mods and candidate["unresolved"] != {}:

        #If mod are still considered but unresolved span is exceeding limit, cut branch
        if candidate["unresolved"]["span"] > mod_coverage_tol and candidate["ac_mod"] <= max_no_mod and candidate["ac_mod_mass"] <= max_mod_mass - min_mod_mass + t_tol:
            new_candidate["Continue"] = False
            #print(new_candidate)
            return new_candidate

        #If matched
        elif match != None:
            #if attempting for mod, add matched mod data to candidate dict and update positions
            if mod != "unmod":
                new_candidate["matches"][candidate["cur_pos"]+1] = [aa, match]
                new_candidate["modifications"][candidate["unresolved"]["pos"]+1] = [candidate["unresolved"]["aa"], mod]
                new_candidate["unresolved"] = {}
                new_candidate["cur_pos"] += 1
                new_candidate["cur_th_mass"] += aa_mass[aa]
                return new_candidate

            #if not attempting for mod, add matched unmod data to candidate dict and update positions
            else:
                new_candidate["matches"][candidate["cur_pos"]+1] = [aa, match]
                new_candidate["cur_pos"] += 1
                new_candidate["cur_th_mass"] += aa_mass[aa]
                new_candidate["unresolved"] = {}
                return new_candidate

        #If not matched to anything, update positions and unresolved info
        else:
            new_candidate["cur_pos"] += 1
            new_candidate["cur_th_mass"] += aa_mass[aa]
            new_candidate["unresolved"]["span"] += 1
            return new_candidate
    
    #If current aa could not be modified
    elif aa not in residue_mods:
        #If matched, add matched unmod data to candidate dict and update positions
        if match != None:
            new_candidate["matches"][candidate["cur_pos"]+1] = [aa, match]
            new_candidate["cur_pos"] += 1
            new_candidate["cur_th_mass"] += aa_mass[aa]
            return new_candidate
        #If not matched, update positions
        else:
            new_candidate["cur_pos"] += 1
            new_candidate["cur_th_mass"] += aa_mass[aa]
            return new_candidate
    
    
#This is a recursive function to perform match and make decision on each position. Checkpoints are set to ensure correct moves
#If multiple modification is matched, mulitple candidate dict branch will be added to a list and they will be searched differently
#When the position reaches the end, the candidate dict will be appended to matched candidate list
def search(candidate, spectrum, sequence, max_mod_mass, basic_mass, fixed_mod_mass, obs_mass, matched_candidates):
    
    #Checkpoint1: check if the sequence end have been reached or branch needed to be cut
    if candidate["cur_pos"] <= len(sequence)-1 and candidate["Continue"] != False:
        #checkpoint passed, get current aa
        aa = sequence[candidate["cur_pos"]]
        
        #Checkpoint 2: If current aa must be modified
        if aa in fixed_mods:
            #Attempt match, make decision and move on
            match = match_ion(17.025965, candidate, aa, fixed_mods[aa], spectrum)
            search(decision(match, candidate, aa, fixed_mods[aa], max_mod_mass), spectrum, sequence, max_mod_mass, basic_mass, fixed_mod_mass, obs_mass, matched_candidates)

        #Checkpoint 3.0: Search n-terminal mods
        elif candidate["cur_pos"] == 0 and len(n_mods) != 0 and candidate["ac_mod"] <= max_no_mod and candidate["ac_mod_mass"] <= max_mod_mass - min_mod_mass + t_tol:
            to_search = []
            #try all mods
            for mod in n_mods:
                #Check if mass after adding the mod will exceed max mass
                if candidate["ac_mod_mass"] + mod_mass[mod] <= max_mod_mass:
                    match = match_ion(17.025965, candidate, aa, mod, spectrum)
                    to_search.append(decision(match, candidate, aa, mod, max_mod_mass))

            #THIS IS THE MAJOR MOD SEARCH NODE. Search for each branch candidate dict.
            for i in to_search:
                search(i, spectrum, sequence, max_mod_mass, basic_mass, fixed_mod_mass, obs_mass, matched_candidates)
        
        #Checkpoint3: If current aa could be modified, and max no of mod has not reached.
        elif aa in residue_mods and candidate["ac_mod"] <= max_no_mod and candidate["ac_mod_mass"] <= max_mod_mass - min_mod_mass + t_tol:
            to_search = []
            #try all mods
            for mod in residue_mods[aa]:
                #Check if mass after adding the mod will exceed max mass
                if candidate["ac_mod_mass"] + mod_mass[mod] <= max_mod_mass:
                    match = match_ion(17.025965, candidate, aa, mod, spectrum)
                    to_search.append(decision(match, candidate, aa, mod, max_mod_mass))

            #THIS IS THE MAJOR MOD SEARCH NODE. Search for each branch candidate dict.
            for i in to_search:
                search(i, spectrum, sequence, max_mod_mass, basic_mass, fixed_mod_mass, obs_mass, matched_candidates)
        
        #Checkpoint4: If aa could not be modified, but there is previously unresolved residue
        elif aa not in residue_mods and candidate["unresolved"] != {}:
            match = match_ion(17.025965, candidate, aa, "unmod", spectrum)
            search(decision(match, candidate, aa, candidate["unresolved"]["mod"], max_mod_mass), spectrum, sequence, max_mod_mass, basic_mass, fixed_mod_mass, obs_mass, matched_candidates)
            
        #Checkpoint5: Perform straight search if aa could not be modified and there is no unresolved mod, or max no or mass of mod has reached
        elif aa not in residue_mods or candidate["ac_mod"] > max_no_mod or candidate["ac_mod_mass"] > max_mod_mass - min_mod_mass:
            match = match_ion(17.025965, candidate, aa, "unmod", spectrum)
            search(decision(match, candidate, aa, "unmod", max_mod_mass), spectrum, sequence, max_mod_mass, basic_mass, fixed_mod_mass, obs_mass, matched_candidates)
        
    #Checkpoint6: Sequence end has reached. Apply post search filtering and append candidate dict to matched candidate list.
    elif candidate["cur_pos"] > len(sequence)-1:
        #print(candidate)
        #Calculate final theoretical mass
        candidate["final_th_mass"] = basic_mass + candidate["ac_mod_mass"] + fixed_mod_mass
        #If final theoretical mass is within tolerance range and can meet min fragment requirement, append to candidate
        if abs(candidate["final_th_mass"] - obs_mass) <= t_tol and len(candidate["matches"]) >= min_frag:
            
            #Add a key for brno proteoform name
            candidate['name'] = ''
            for mod in candidate["modifications"]:
                candidate["name"] += candidate["modifications"][mod][0] + str(mod) + '(' + candidate["modifications"][mod][1] + ')' + " "
            #Add a key for modification state
            #candidate["mod_state"] = (candidate["final_th_mass"], candidate["name"].count("ac"))
            matched_candidates.append(candidate)