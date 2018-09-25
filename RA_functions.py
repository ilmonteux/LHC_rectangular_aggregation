import time, sys, os
from joblib import Parallel, delayed # for parallel processing
import numpy as np
import math
from itertools import combinations
from interval import Interval as interval

N_PARALLEL_JOBS = 4 # it saturates the CPU with N=4 on my mac

# variables that will require closed intervals
int_variables = []

# reorder kinematic variables to have Nx variables appear first, MET/MT2 last
def reorder_vars(kv):
    ''' Defines 'nice' order of kinematic variables, e.g. Nj, Nb, HT, MET, ... At same time, save list of integer variables'''
    if kv.startswith('N'): 
        pos = 10
        if kv == 'Nj': pos = 1
        if kv not in int_variables: int_variables.append(kv)
    # add other ifs as needed
    
    # put MET, MT2 as last
    elif kv in ['MET','MHT','MT2']: 
        pos = 100
    # everything else has same priority unless specified above
    else: 
        pos = 50
    return pos

# deals with overflow bins, where the max=Inf.
def getn(num):
    '''returns numerical value from numeric string. If string = Inf. (one-sided bin), returns large number (99.999 TeV)'''
    if num != 'Inf.': return float(num)
    else: return 99999. # overflow

def read_SR_defs(input_file):
    ''' Reads signal region definitions from a space- or tab-delimited file, in the following format:
      Bin\t KV_1 \t\t KV_1 \t\t KV_2 \t\t KV_2 \t\t...KV_n... \t (Bkg.\t ...)
      bin\t KV_1^min\t KV_1^max\t KV_2^min\t KV_2^max\t...n... \t(Bkg.\t...(bkg errors)...\tData)
    '''
    file = open(input_file).readlines()
    # opens data file, reads
    kin_vars=[]
    kin_map = {}
    err_types=[]
    err_map={}
    # variables are listed in header of the file
    for i_col, var in enumerate(file[0].split()):
        # do not read event counts
        if var == 'Bkg': bg_column = i_col
        elif var == 'Data': obs_column = i_col
        elif var.startswith('bg_err'):
            err_types.append(var)
            err_map[var] = i_col
        # skip bin number, assuming it's first column
        elif var != 'Bin':
            # if repeated kin. variable, two columns are assumed to be min/max for that variable
            if len(kin_vars)>0 and var in kin_vars:
                kin_map[var]=[kin_map[var][0],i_col]
            else:
                kin_map[var] = [i_col]
                kin_vars.append(var)

    # bg errors: assume they are +stat+sys,-stat-sys
    if len(err_types) == 4: err_func = lambda x: np.sqrt(max([ x[0]**2 + x[1]**2, x[2]**2 + x[3]**2 ]))
    # bg errors: assume they are +stat,-stat,sys
    if len(err_types) == 3: err_func = lambda x: np.sqrt(max([ x[0]**2 + x[2]**2, x[1]**2 + x[2]**2 ]))
    # bg errors: assume they are +err,-err
    if len(err_types) == 2: err_func = lambda x: np.sqrt(max([ x[0]**2, x[1]**2 ]))
    # bg error
    if len(err_types) == 1: err_func = lambda x: x[0]

    # redefine sorted kinematic variables
    kin_vars = sorted(kin_vars, key=reorder_vars)

    print '\nRunning rectangular aggregation algorithm on given input...'
    print 'Recognized kinematic variables:', kin_vars
    print 'Mapping between variable and column in data file:', kin_map
    
    # defined global so they do not need to be passed as arguments to all functions
    global var_min, var_max 
    var_min, var_max = {}, {}
    binlist = []
    obs_map, bg_map, bg_err_map = {}, {}, {}
    # now read following lines to get min/max of each bin in each variable
    for line in file[1:]:
        ll = line.split()
        # read bin number and add to binlist
        Bin = int(ll[0])
        binlist.append(Bin)
        # read rest of line and defines min/max mapping from each column
        var_min[Bin]={}
        var_max[Bin]={}
        for var in kin_vars:
            # read column element corresponding to each KV. If single element, doubles it as min=max
            var_delimiter = [ll[i] for i in kin_map[var]]
            if len(var_delimiter)==1: var_delimiter = [var_delimiter[0], var_delimiter[0]]
            var_min[Bin][var] = var_delimiter[0]
            var_max[Bin][var] = var_delimiter[1]

        # get observed counts, backgrounds and errors
        bg_err_map[Bin] = err_func([ float(ll[err_map[e]]) for e in err_types ])
        obs_map[Bin] = int(ll[obs_column])
        bg_map[Bin] = float(ll[bg_column])
    
    # return event counts taking care of nasty case in which bins are not defined sequentially in the file
    binlist = sorted(binlist)
    obs_vec = [obs_map[Bin] for Bin in sorted(obs_map)]
    bg_vec = [bg_map[Bin] for Bin in sorted(bg_map)]
    bg_err_vec = [bg_err_map[Bin] for Bin in sorted(bg_err_map)]

    return [binlist, obs_vec, bg_vec, bg_err_vec, kin_vars, kin_map, var_min, var_max]

# function asking if a certain bin is in a range for a kinematic variable
def isbininrange(Bin, var, range):
    '''Is the bin overlapping with a certain range in a given variable, given bin definitions var_min/max?'''
    # if var=Nx, it is integer and the interval range is a closed interval, otherwise it's open
    closed_interval = True if var.startswith('N') else False
    # define interval
    range_interval = interval(getn(range[0]),getn(range[1]),lower_closed=closed_interval, upper_closed=closed_interval)
    
    # does the bin overlap with the interval?
    res = ( (getn(var_min[Bin][var]) in range_interval) or (getn(var_max[Bin][var]) in range_interval) or (getn(var_min[Bin][var]) <= getn(range[0]) and getn(var_max[Bin][var])>=getn(range[1])))

    return res

# function that makes an RA given a list of variables and ranges for each variable
def makeRA(var_range_list):
    '''Select bins inside an aggregation defined by a list of kinematic variables and corresponding ranges that define the RA'''
    vars = [var_range[0] for var_range in var_range_list]
    ranges = [var_range[1] for var_range in var_range_list]
    binlist = var_min.keys()
    return [ Bin for Bin in binlist if np.prod([isbininrange(Bin,var,range) for var,range in zip(vars, ranges)]) ]

# function to print bins included in RA
def printRA(RA,dictionary):
    ''' Print definition of each bin given a certain RA (or list of bin numbers)'''
    for Bin in RA:
        i_bin = [ dic[0][1] for dic in dictionary].index(Bin)
        print dictionary[i_bin]

def get_RA_from_bins(binlist, dictionary):
    """Returns range in each kinematic variable that is covered by the input list of bins"""
    mins = {}
    maxs = {}
    ranges = {}
    kin_vars = sorted(dictionary.values()[0].keys(), key=reorder_vars)
    for var in kin_vars:
        m = min([dictionary[Bin][var][0] for Bin in binlist ])
        M = max([dictionary[Bin][var][1] for Bin in binlist ])
        ranges[var] = (m,M)
    return ranges
    
def varbins(var):
    ''' For a given kinematic variable, return all possible 1D intervals between SR separators'''
    binlist = var_min.keys()
    delimiters = sorted(list(set([ (var_min[Bin][var]) for Bin in binlist])),key=getn)
    delimiters2 = sorted(list(set([ (var_max[Bin][var]) for Bin in binlist])),key=getn)
    if var in int_variables:
        res = [ [n1,n2] for n1 in delimiters for n2 in delimiters2 if getn(n1)<= getn(n2)]
    else:
        delimiters = sorted(list(set([ (var_min[Bin][var]) for Bin in binlist] + [ (var_max[Bin][var]) for Bin in binlist])),key=getn)
        res = [ [n1,n2] for n1 in delimiters for n2 in delimiters if getn(n1)< getn(n2)]
    return res

# routine defining a rectangular aggregation:
# recursive because length of kinematic variables depends on input
def recursive_RA_loop(RAdef, vars):
    ''' Recursive method to find all possible rectangular aggregations given a list of kinematic variables'''
    var = vars[0]
    RAdef0 = RAdef
    res = []
    if len(vars)==1:
        for var_range in varbins(var):
            res.append(RAdef0 + [[var,var_range]])
        return res
    else:
        for var_range in varbins(var):
            RAdef = recursive_RA_loop(RAdef0 + [[var,var_range]],vars[1:])
            res=res+(RAdef)
        return res
    
# simple seconds > hour_minute_seconds
def sec_to_hms(t):
    t=int(t)
    out=""
    h, m, s = t//3600, (t%3600)//60, (t%3600)%60
    if h > 0: out=out+str(h)+"h"
    if m > 0: out=out+str(m)+"m"
    out=out+str(s)+"s"
    return out


# grouping function that collects overlapping aggregations
def grouper(sequence):
    result = []  # will hold (members, group) tuples

    for item in sequence:
        for members, group in result:
            # if members.intersection(item):  # any overlap within lists
            # lists overlap by at least half of elements
            if len(members.intersection(item)) >= 0.5*min(map(len,[members,item])): 
                members.update(item)
                group.append(item)
                break
        else:  # no group found, add new
            result.append((set(item), [item]))
    return [group for members, group in result]