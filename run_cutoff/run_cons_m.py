import numpy as np

mass=40;ncore=32
arrlen=100
ui=-2;uf=2
q=1;

eref1=np.array([1.1,1.15,1.2,1.25])
bref1=np.array([50,60,80,100]) 

# eref1=np.array([1.1,1.15,])
# bref1=np.array([50,60]) 

ecc1=np.linspace(1.09,1.4,arrlen)

impacts1=np.linspace(40,150,arrlen)

import pycbc.waveform
from pycbc.waveform import td_approximants
import os
import matplotlib.pyplot as plt
import pylab
from pycbc.filter import match
from pycbc.psd import aLIGOZeroDetHighPower
#from constants import *
from pycbc.waveform import get_td_waveform
from multiprocess import Pool
import pickle


import os
fn='file_m_'+str(mass)+'_hc2_LAL.pkl'
os.remove(fn) if os.path.exists(fn) else None  #here fn - filename or filepath

#eref,bref,eccs,impacts=np.meshgrid(eref1,bref1,ecc1,impacts1)

approximant = 'HyperbolicTD'


delta_t = 1./1024
inc=np.pi/3
#distance=1e6*pc
order=3
f_lower = 10
matchfactor=0.99

def cal_match_eb(eref1,bref1,e1,b1):
    
    #hp,hc=get_hyp_waveform(mref1,q,ecc, bref1,  delta_t, inc, distance, order)

    hp, hc = get_td_waveform(approximant='HyperbolicTD',
                                 mass1=mass/2,
                                 mass2=mass/2,
                                 delta_t=delta_t,
                                 hyp_eccentricity = eref1,
                                 b = bref1,
                                 inclination = inc,
                                 ui = ui,
                                 uf = uf,
                                 distance = 1,
                                 phase_order = -1,
                                 f_lower = f_lower)

    psd = aLIGOZeroDetHighPower(len(hc) // 2 + 1, 1.0 / hc.duration, f_lower)
    
    hp2, hc2 = get_td_waveform(approximant='HyperbolicTD',
                                 mass1=mass/2,
                                 mass2=mass/2,
                                 delta_t=delta_t,
                                 hyp_eccentricity = e1,
                                 b = b1,
                                 inclination = inc,
                                 ui = ui,
                                 uf = uf,
                                 distance = 1,
                                 phase_order = -1,
                                 f_lower = f_lower)
    hc2 = hc2[:len(hc)] if len(hc) < len(hc2) else hc2
    hc2.resize(len(hc))

    m, idx = match(hc, hc2, psd=psd, low_frequency_cutoff=f_lower)
    return m

# arr2=np.zeros((len(mref1),len(bref1),arrlen,arrlen))

peakfile = open(fn, 'ab')

from multiprocess import Pool

if __name__ == "__main__":

    pool = Pool(ncore)
    results = pool.map(lambda t : cal_match_eb(t[0],t[1],t[2],t[3]), 
                       [(eref1[i],bref1[j],ecc1[k],impacts1[l])
                        for i in range(len(eref1)) for j in range(len(bref1)) 
                        for k in range(arrlen) for l in range(arrlen) ] )
    pool.close()
    pickle.dump(np.array(results).reshape((len(eref1),len(bref1),arrlen,arrlen)),peakfile)
    peakfile.close()



