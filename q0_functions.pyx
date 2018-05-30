# compile this file with:
# python cython_setup.py build_ext --inplace

cimport cython
from iminuit import Minuit, describe, Struct
import numpy as np
cimport numpy as np
from scipy.linalg import block_diag
from libc.math cimport log, M_PI
from cpython cimport array
from libcpp cimport bool
import array
from collections import OrderedDict



cdef array.array obs_vec
cdef array.array back_vec
cdef array.array backerr_vec
cdef np.ndarray Vcovar
cdef np.ndarray inverseV

#function to load event counts
def loadObs(obs):
    global obs_vec0
    obs_vec0 = array.array('d', obs )

# function to load background (expected) vector
def loadBack(back):
    global back_vec0
    back_vec0 = array.array('d', back )

# function to load background error vector
def loadBackErr(back_err):
    global backerr_vec0
    backerr_vec0 = array.array('d', back_err )

# function to load all event counts
def loadAllCounts(obs, back, back_err):
    global obs_vec0, back_vec0, backerr_vec0
    obs_vec0 = array.array('d', obs )
    back_vec0 = array.array('d', back )
    backerr_vec0 = array.array('d', back_err )
    return obs_vec0, back_vec0, backerr_vec0

# function to load covariance matrix
@cython.embedsignature(True)
def loadVcovar(path_to_V):
    global Vcovar
    try: Vcovar = np.loadtxt(path_to_V)
    except: 
        Vcovar = np.diag([err**2 for err in backerr_vec0 ])
        print "No covariance matrix: will build diagonal V from single bin uncertaintites..."

def setAll():
    global obs_vec, back_vec, backerr_vec, Vcovar, inverseV
    obs_vec, back_vec, backerr_vec = loadAllCounts(obs_vec0, back_vec0, backerr_vec0)
    inverseV = np.linalg.inv(Vcovar)
    
# function to load signal vector
cdef array.array sig_vec
cpdef np.int_t setSig(np.ndarray[np.double_t, ndim=1] sig):
    global sig_vec
    sig_vec = array.array('d', sig )
    return 1

# SR aggregation code: input is list of bins (in mathematica format, 1=first bin, will get translated to python->0)
def aggregate(SR_list):
    # load search counts, before merging bins into aggregation
    global obs_vec, back_vec, backerr_vec, Vcovar, inverseV
    obs_vec, back_vec, backerr_vec = loadAllCounts(obs_vec0, back_vec0, backerr_vec0)
    
    # translate from 1-based human list to 0-based python list
    SR_list = [SR-1 for SR in SR_list]
    out_SR_list = [ n for n in range(len(obs_vec)) if n not in SR_list ]

    # create reduced vectors of obs, backgrounds, errors
    obs_ASR = sum([obs_vec[SR] for SR in SR_list])
    obs_vec = array.array('d',[obs_ASR] + [obs_vec[SR] for SR in out_SR_list])
    back_ASR = sum([back_vec[SR] for SR in SR_list])
    back_vec = array.array('d',[back_ASR] + [back_vec[SR] for SR in out_SR_list])
    backerr_ASR = np.sqrt(sum([backerr_vec[SR]**2 for SR in SR_list]))
    backerr_vec = array.array('d',[backerr_ASR] + [backerr_vec[SR] for SR in out_SR_list])
    # create new covariance matrix
    V11 = sum(Vcovar[SR1][SR2] for SR1 in SR_list for SR2 in SR_list )
    V1i = [sum(Vcovar[SR1][SR2] for SR1 in SR_list )  for SR2 in out_SR_list]
    Vij = [[Vcovar[SR1][SR2] for SR1 in out_SR_list] for SR2 in out_SR_list ]
    # protect against extreme scenario where all SRs are one aggregation
    if len(V1i)>0:
        newV = np.bmat([[np.matrix(V11),np.matrix(V1i)],[np.transpose(np.matrix(V1i)),np.matrix(Vij)]])
    else: newV=[[V11]]
    inverseV = np.linalg.inv(newV)
    
    return [obs_vec, back_vec, newV]


def call_Minuit(function, **kwargs):
    # setting up default arguments for calling minuit: 
    #   - all theta_i's have lower limit -back[i] 
    #   - signal events are input as fixed

    # variables that will be minimized over: signal multiplier mu, nuisance parameters theta_i
    param_names = ['mu'] + [ 't'+str(t) for t in range(len(back_vec))]
    #t_names = [ 't'+str(t) for t in range(len(back_vec))]
    
    ## initialize nuisance parameters at zero
    #th_in = [0. for bg in back_vec]
    # nuisance parameters are limited by bg+th > 0 (otherwise the Poisson statistics doesnot make sense)
    # restrict minimization to only the physical region
    th_limit_names = [ 'limit_t'+str(t) for t in range(len(back_vec))]
    th_limits = [(-bg,None) for bg in back_vec]
    args_default = [('forced_parameters',param_names),('print_level',-1),('pedantic',False)]+zip(th_limit_names, th_limits)

    # read additional kwargs for minuit call - append them to the variables just defined
    new_dict = OrderedDict(args_default)
    for key, value in kwargs.items():
        new_dict[key] = value
    return Minuit(fcn=function, **new_dict)




#this is slightly faster than math.factorial
cdef np.double_t myfac(np.double_t n):
    if n>1: return myfac(n-1.)*n
    else: return 1.

cdef np.double_t logfac(np.double_t n): # Ramanujan's formula
    if n > 10:  return n*log(n) - n + (log(n*(1.+4.*n*(1.+2.*n))))/6. + log(M_PI)/2.
    else: return log(myfac(n))

#evaluate LL at given mu and theta vector
# LL = -2 log( Poisson(n_obs, mu*signal + bg + theta) exp(-1/2 theta V^-1 theta)
@cython.boundscheck(False)
cdef np.double_t LLfunc1(np.double_t mu, np.ndarray[np.double_t, ndim=1] theta_vec):
    cdef np.double_t ll = -0.5 * np.dot(np.dot(theta_vec,inverseV),theta_vec)
    cdef np.double_t lam
    cdef Py_ssize_t i
    for i in range(0,len(obs_vec)):
        lam = mu * sig_vec[i] + back_vec[i] + theta_vec[i]
        if lam<=0: lam = 10**(-8.)
        ll+= obs_vec[i] * log(lam) - lam - logfac(obs_vec[i])
    return -2.*ll

#rearranges inputs of LLfunc1 into form which Minuit can call
def LL1(mu,*tvec):
    return LLfunc1(mu, np.array(tvec,np.double))


#runs minuit to find minimum LL
# was cpdef (np.double_t,np.double_t) and it was crashing
@cython.embedsignature(True)
cpdef list runmin(np.double_t tmu, bint flag_fix_mu):
    mLL= call_Minuit(LL1, mu=tmu, fix_mu = flag_fix_mu, print_level = 0)

    cdef bint flag_done = False
    cdef np.double_t edm = 1.
    cdef np.double_t func_val = 0.
    while not flag_done:
        mLL.migrad()
        flag_done = mLL.get_fmin().is_valid
        if abs(mLL.edm/edm-1) < 0.001: flag_done = True
        edm = mLL.edm
        if abs(func_val/mLL.fval -1)<.0001: flag_done = True
        func_val = mLL.fval 
        
    
    # returns parameter values (mu, theta_i) at minimum, and value of function
    return [ mLL.args[0], mLL.fval ]

#same as above, but also return post-fit error(s) after minimization
cdef array.array postfit_err_vec
@cython.embedsignature(True)
cpdef list runmin_err(np.double_t tmu, bint flag_fix_mu):
    mLL= call_Minuit(LL1, mu=tmu, fix_mu = flag_fix_mu, print_level = 0)

    cdef bint flag_done = False
    cdef np.double_t edm = 1.
    cdef np.double_t func_val = 0.

    # keep minimizing until minimum is found, which is when:
    #   * Minuit says it has found the minimum (get_fmin().is_valid)
    #   * or, when successive evaluations do not change the minimzation:
    #       - when the log-likelihood value after each evaluation stays the same (at the 10^-4 level)
    #       - when edm (Estimated distance to minimum) is not changing
    while not flag_done:
        mLL.migrad() # minimize
        
        flag_done = mLL.get_fmin().is_valid
        if abs(mLL.edm/edm-1) < 0.001: flag_done = True
        edm = mLL.edm
        if abs(func_val/mLL.fval -1) < .0001: flag_done = True
        func_val = mLL.fval 
    
    t_names = [ 't'+str(t) for t in range(len(back_vec))] # labels of theta vector
    postfit_err_vec = array.array('d', [ mLL.errors[th] for th in t_names ] )

    # returns parameter values (mu, theta_i) at minimum, value of function, postfit theta's and relative errors
    return [mLL.args[0], mLL.fval, mLL.args[1:], postfit_err_vec ]

