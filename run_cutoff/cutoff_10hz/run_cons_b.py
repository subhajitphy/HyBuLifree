import numpy as np

impact=60;ncore=32
arrlen=100
ui=-2;uf=2;
q=1;

mref1=np.array([20,40,60,100])
eref1=np.array([1.1,1.15,1.2,1.25])

masses1=np.linspace(10,150,arrlen)

ecc1=np.linspace(1.09,1.4,arrlen)

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
fn='file_b_'+str(impact)+'_hc2_LAL.pkl'
os.remove(fn) if os.path.exists(fn) else None  #here fn - filename or filepath

#mref,eref,masses,eccs=np.meshgrid(mref1,eref1,masses1,ecc1)

approximant = 'HyperbolicTD'


delta_t = 1./1024
inc=np.pi/3
#distance=1e6*pc
order=3
f_lower = 10
matchfactor=0.99

def cal_match_me(mref1,eref1,M1,e1):
    
    #hp,hc=get_hyp_waveform(mref1,q,ecc, bref1,  delta_t, inc, distance, order)

    hp, hc = get_td_waveform(approximant='HyperbolicTD',
                                 mass1=mref1/2,
                                 mass2=mref1/2,
                                 delta_t=delta_t,
                                 hyp_eccentricity = eref1,
                                 b = impact,
                                 inclination = inc,
                                 ui = ui,
                                 uf = uf,
                                 distance = 1,
                                 phase_order = -1,
                                 f_lower = f_lower)

    psd = aLIGOZeroDetHighPower(len(hc) // 2 + 1, 1.0 / hc.duration, f_lower)
    
    hp2, hc2 = get_td_waveform(approximant='HyperbolicTD',
                                 mass1=M1/2,
                                 mass2=M1/2,
                                 delta_t=delta_t,
                                 hyp_eccentricity = e1,
                                 b = impact,
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
    results = pool.map(lambda t : cal_match_me(t[0],t[1],t[2],t[3]), 
                       [(mref1[i],eref1[j],masses1[k],ecc1[l])
                        for i in range(len(mref1)) for j in range(len(eref1)) 
                        for k in range(arrlen) for l in range(arrlen) ] )
    pool.close()

    pickle.dump(np.array(results).reshape((len(mref1),len(eref1),arrlen,arrlen)),peakfile)
    peakfile.close()


