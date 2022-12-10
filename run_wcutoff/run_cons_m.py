import numpy as np

mass=40;ncore=32
arrlen=100
ui=-2;uf=2
q=1;

eref1=np.array([1.1,1.15,1.2,1.25])
bref1=np.array([50,60,80,100]) 


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

from hyp_td_waveform import get_hyp_waveform
from constants import *

fn='file_m_'+str(mass)+'_hc2_LAL.pkl'
os.remove(fn) if os.path.exists(fn) else None  #here fn - filename or filepath

#eref,bref,eccs,impacts=np.meshgrid(eref1,bref1,ecc1,impacts1)

approximant = 'HyperbolicTD'


delta_t = 1./1024
inc=np.pi/3
#distance=1e6*pc
order=3

matchfactor=0.99

def cal_match_eb(eref1,bref1,e1,b1):
    

    hp,hc,uarr=get_hyp_waveform(mass,q,eref1,bref1,delta_t,inc,1e6*pc,order,"dudt")

    psd = aLIGOZeroDetHighPower(len(hc) // 2 + 1, 1.0 / hc.duration, 10)
    
    
    hp2,hc2,uarr2=get_hyp_waveform(mass,q,e1,b1,delta_t,inc,1e6*pc,order,"dudt")

    hc2 = hc2[:len(hc)] if len(hc) < len(hc2) else hc2
    hc2.resize(len(hc))

    m, idx = match(hc, hc2, psd=psd, low_frequency_cutoff=10)
    return m

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




