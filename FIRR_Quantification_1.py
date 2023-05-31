# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 15:27:56 2021

Functions to check PrSMs SDIs and perform co-isolated PrSM quantification

@author: wingkinlui
"""

import statistics
from scipy import stats
from pulp import *

"""
Parameters
"""
nIB = [["ac","me,me,me"]]
mod_coverage_tol = 5
noise = 0.01

#A holder to contain PulP variables
variables_holder = ["c" + str(i) for i in range(0, 1000)]


def informative_partition(current_group, mutable_sites, informative_groups, ):
#The recursive function to partition informative groups by SDI coincidences

    #Only progress when the site is before or is the second last site
    if current_group["cur_site"] < mutable_sites[-1]:
        #print("Current site: ")
        #print(current_group["cur_site"])
        
        #Get all possible SDIs
        SDIs = [i for i in range(current_group["cur_site"], mutable_sites[mutable_sites.index(current_group["cur_site"])+1])]

        #Get total no. of remained candidates in group
        n_candidates = len(current_group["candidates"])

        #Try all SDI
        for SDI in SDIs:
            coincidence = []
            #Get candidates having the SDI and calculdate SDI conincidence rate
            for candidate_id in current_group["candidates"]:
                if SDI in list(candidates[candidate_id]["matches"].keys()):
                    coincidence.append(candidate_id)
            coincidence_rate = len(coincidence)/n_candidates

            
            #If more than one candidate has such SDI, a new branch of informative group is created
            if coincidence_rate > 1/n_candidates:
                new_group = {}

                #Update mutable site of the new group
                mutable_sites = []
                for candidate_i in coincidence:
                    mutable_sites += candidates[candidate_i]["modifications"].keys()
                mutable_sites = sorted(list(set(mutable_sites)))
                
                #Update site
                #If the current site is already later than the last site of the new group, then use the current site and let it enter else condition in the next recursion
                if current_group["cur_site"] >= mutable_sites[-1]:
                    new_cur_site = current_group["cur_site"]
                else:
                    #Otherwise, find the next site relative to the new group
                    for site in mutable_sites:
                        if site > current_group["cur_site"]:
                            new_cur_site = site
                            break

                new_group["cur_site"] = new_cur_site
                new_group["candidates"] = coincidence
                informative_partition(new_group, mutable_sites, informative_groups)

    #When the site reaches the last, need to check if isosite near-isobaric mod is assigned there
    else:
        #Get what the last mutable position modifications are
        last_mod = []
        for candidate_i in current_group["candidates"]:
            if mutable_sites[-1] in candidates[candidate_i]["modifications"].keys():
                last_mod.append(candidates[candidate_i]["modifications"][mutable_sites[-1]][1])

        #Check if near-isobaric mod is present
        for nIB_set in nIB:
            if len(set(nIB_set).intersection(set(last_mod))) > 1:
                #All possible ions between the last site and a coverage tolerance range beyond will be tried
                SDIs = [i for i in range(mutable_sites[-1], mutable_sites[-1] + mod_coverage_tol+1)]
                n_candidates = len(current_group["candidates"])

                for SDI in SDIs:
                    coincidence = []

                    #Get candidates having the SDI and calculdate SDI conincidence rate
                    for candidate_id in current_group["candidates"]:
                        if SDI in list(candidates[candidate_id]["matches"].keys()):
                            coincidence.append(candidate_id)

                    coincidence_rate = len(coincidence)/n_candidates

                    #If more than one candidate has post-end site SDI, the informative groups is added to the result list
                    if coincidence_rate > 1/n_candidates:
                        informative_groups.append(tuple(coincidence))

                #End site near-isobaric modification will only be checked once
                break

            else:
                #If there is no end site near-isobaric mod, just add the resulted informative group into result list
                informative_groups.append(tuple(current_group["candidates"]))


def FIRR_dict_create(candidates, scan, report):
#This is a function to create a FIRR dictionary from list of candidates
    #Update mutable sites in case some is filtered in partitioning
    mutable_sites = []
    for candidate in candidates:
        mutable_sites += candidate["modifications"].keys()
    mutable_sites = sorted(list(set(mutable_sites)))
    
    FIRR_dict = {}
    for site_i in range(len(mutable_sites)):
        acc_mod_states = {}

        #Get acc_mod_state and respective candidate compositions
        for candidate_i, candidate in enumerate(candidates):
            mod_state = [candidate["modifications"][pos][1] for pos in candidate["modifications"].keys() if pos <= mutable_sites[site_i]]
            if len(mod_state) < 1:
                mod_state = ["unmod"]
            mod_state = ",".join(sorted(mod_state))
            if mod_state not in acc_mod_states:
                acc_mod_states[mod_state] = [[candidate_i],[]]
            else:
                acc_mod_states[mod_state][0].append(candidate_i)

        #Site with all candidates having the same mod state will be ignored;
        #This applies to state in the middle where all candidate temporary have the same state
        #This also applies to the end of the mutable site where all candidate have the same state;
        #Except if there are isosite near-isobaric modification
        if len(acc_mod_states[mod_state][0]) == len(candidates):
            continue

        #Get common SDIs
        SDIs = []
        for candidate in candidates:
            #In normal case, SDI in-between mutable sites will be investigated
            if site_i < len(mutable_sites)-1:
                SDIs.append([i for i in candidate["matches"].keys() if (i >= mutable_sites[site_i] and i < mutable_sites[site_i+1])])

            #In the case where there are isosite near-isobaric modification, the end-site mod state will still not be the same. Ion following the site will be considered
            else:
                SDIs.append([i for i in candidate["matches"].keys() if i >= mutable_sites[site_i]])
        
        common_SDIs = set.intersection(*map(set,SDIs))

        #If the isosite near-isobaric modification does not occur at the end-position, post-end position SDI are dispensible.
        #In such case, if there are no post-end position SDI, its OK to ignore them.
        if len(common_SDIs) == 0:
            continue
        
        #Calculate FIRR:
        for SDI in common_SDIs:
            
            fragments = []
            for state in acc_mod_states:
                fragments.append(candidates[acc_mod_states[state][0][0]]["matches"][SDI][1])
            
            #Isosite near-isobaric modification may lead to converged ions, leading to identical SDI. In that case, ignore the SDI. Subsequent SDI would also be ignored due to convergence.
            if len(set(fragments)) != len(fragments):
                continue
            
            #Calculate FIRR for each position
            total_FI = sum([i[1] for i in fragments])
            for fragment, state in zip(fragments, acc_mod_states):
                acc_mod_states[state][1].append(fragment[1]/total_FI)
        
        #If all the SDI in the site are converged ions, the FIRR list will be empty. Prevent the site to be appended to the FIRR dict.
        if len(acc_mod_states[state][1]) == 0:
            report["Converged"].append(scan)
            continue
        
        #Drop SDI with FIRR > mean +/- deviation
        for state in acc_mod_states:
            if len(acc_mod_states[state][1]) > 1:
                mean = statistics.fmean(acc_mod_states[state][1])
                sd = statistics.stdev(acc_mod_states[state][1])
                acc_mod_states[state][1] = [i for i in acc_mod_states[state][1] if (i >= mean-sd and i<= mean+sd)]
        
        FIRR_dict[mutable_sites[site_i]] = acc_mod_states
    return FIRR_dict

def FIRR_flow_processing(FIRR_dict, candidates):
    #This is the FIRR flow filtering function. It takes in candidates and corresponding FIRR dict. Returns filtered candidates and FIRR dict
    #Update mutable sites
    mutable_sites = [i for i in FIRR_dict.keys()]
    
    CIMs = []
    for site_i in range(len(mutable_sites[:-1])):

        #Get mod that are present in both consecutive site
        unchanged_mod = set(FIRR_dict[mutable_sites[site_i]].keys()).intersection(set(FIRR_dict[mutable_sites[site_i+1]].keys()))

        #Observe FIRR change and get putative CIM
        for mod in unchanged_mod:
            pre = FIRR_dict[mutable_sites[site_i]][mod][1]
            post = FIRR_dict[mutable_sites[site_i+1]][mod][1]
            pre_compo = set(FIRR_dict[mutable_sites[site_i]][mod][0])
            post_compo = set(FIRR_dict[mutable_sites[site_i+1]][mod][0])

            #If both modstate has >1 SDIs, use t-test to determine change in FIRR
            if len(pre) > 1 and len(post) > 1:
                if stats.ttest_ind(pre,post,equal_var=False)[1] > 0.05:
                    #print(mutable_sites[site_i])
                    #print(unchanged_mod)
                    CIMs += list(pre_compo.symmetric_difference(post_compo))

            #If one of the modstate has only 1 SDI, assume no change if FIRR change is < defined noise
            else:
                if abs(statistics.fmean(pre) - statistics.fmean(post)) < noise:
                    CIMs += list(pre_compo.symmetric_difference(post_compo))
    
    #Aftering getting the list of CIMs, filter them out from the candidates
    candidates = [candidates[i] for i in range(len(candidates)) if i not in CIMs]
    return candidates

def LP_quant(FIRR_dict, candidates):
    #Holds all coefficient: {site_mod: [candidate composition], mean, deviation}
    site_mod_FIRR = {}
    for site in FIRR_dict:
        for mod in FIRR_dict[site]:
            firr = FIRR_dict[site][mod][1]
            if len(firr) > 1:
                site_mod_FIRR["".join([str(site),mod])] = [FIRR_dict[site][mod][0], statistics.fmean(firr), statistics.stdev(firr)]
            else:
                site_mod_FIRR["".join([str(site),mod])] = [FIRR_dict[site][mod][0], firr[0] , noise]
    
    #Initialize Linear programming model
    LP = LpProblem("Quantification", LpMinimize)
    
    #Set variables
    variables = [LpVariable(variables_holder[candidate_id], lowBound = 0, upBound = 1, cat='Continuous') for candidate_id in range(len(candidates))]
    
    #Objective function
    obj_func = lpSum([(site_mod_FIRR[site_mod][1]-lpSum([variables[var] for var in site_mod_FIRR[site_mod][0]]))/(site_mod_FIRR[site_mod][2]) for site_mod in site_mod_FIRR])
    LP += obj_func
    
    for site_mod in site_mod_FIRR:
        LP += (site_mod_FIRR[site_mod][1]-lpSum([variables[var] for var in site_mod_FIRR[site_mod][0]]) <= site_mod_FIRR[site_mod][2])
        LP += (site_mod_FIRR[site_mod][1]-lpSum([variables[var] for var in site_mod_FIRR[site_mod][0]]) >= -site_mod_FIRR[site_mod][2])
    
    #Solve
    LP.solve()
    
    quant_results = []
    if LpStatus[LP.status] == "Optimal":
        for var in range(len(variables)):
            if variables[var].varValue >= 0:
                candidates[var]["Abundance"] = precursor_ab * variables[var].varValue
                quant_results.append(candidates[var])
        return quant_results

def Quantification(search_results, report):

    #This is the main driver function of post search processing and quantification
    post_process_results = {}

    #Iterate scans in results
    progress = 0
    total_scans = len(search_results)
    for scan in search_results:
        progress += 1
        print("Quantifying " + str(scan)+ "; " + str(progress) + " of " + str(total_scans) + " completed")
        global precursor_ab
        precursor_ab = search_results[scan]["precursor_intensity"]
        informative_groups = []
        #Declare candidates and mutable sites as global for ease of use
        global candidates
        candidates = search_results[scan]["candidates"]
        #No need of quantification if there is only 1 candidate:
        if len(candidates) == 1:
            post_process_results[scan] = search_results[scan]
            post_process_results[scan]["candidates"][0]["Abundance"] = precursor_ab
            continue
        
        #Initiate mutable sites
        mutable_sites = []
        for candidate in candidates:
            mutable_sites += candidate["modifications"].keys()
        mutable_sites = sorted(list(set(mutable_sites)))
        
        #initiate informative groups, first group contains no discard
        current_group = {
            "cur_site": mutable_sites[0],
            "candidates": [i for i in range(len(candidates))]}
        
        
        #Conduct informative group partition
        informative_partition(current_group, mutable_sites, informative_groups, )
        
        
        #Pass scan if there is no informative group
        if len(informative_groups) == 0:
            report["NoInformativeGrp"].append(scan)
            continue
        
        #Get the biggest group
        informative_groups = set(informative_groups)
        max_candidates_no = len(sorted(informative_groups, key = len, reverse= True)[0])
        max_informative_group = [i for i in informative_groups if len(i) == max_candidates_no]
        
        
        #If there is more than one largest group, there is no choice but to discard spectrum
        if len(max_informative_group) > 1:
            continue
        else:
            candidates = [candidates[i] for i in max_informative_group[0]]
        
        #Using the shortlisted informative group, create FIRR dictionary
        FIRR_dict = FIRR_dict_create(candidates, scan, report)
        if len(FIRR_dict) == 0:
            report["All_Filtered"].append(scan)
            continue
            
        '''
        #Conduct FIRR flow investigation and recreate FIRR_dict
        candidates = FIRR_flow_processing(FIRR_dict, candidates)
        FIRR_dict = FIRR_dict_create(candidates, scan, report)
        if len(FIRR_dict) == 0:
            report["All_Filtered"].append(scan)
            continue
       '''
        
        #Conduct Linear programing quantification
        candidates = LP_quant(FIRR_dict, candidates)
        
        if candidates is not None:
            post_process_results[scan] = search_results[scan]
            post_process_results[scan]["candidates"] = candidates
            report["Quantified"].append(scan)
        else:
            report["UndefinedSol"].append(scan)
            continue
            
    return post_process_results, report
        
        