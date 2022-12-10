import matplotlib.pyplot as plt
import pylab
from pycbc.filter import match
import numpy as np

from hyp_td_waveform_v2 import get_hyp_waveform
from constants import *


M=200;q=1
f_lower = 10
delta_t = 1./4096
ecc=1.1
impact=50
inc=pi/3
distance=1e6*pc
order=3

hp,hc=get_hyp_waveform(M,q,ecc, impact,  delta_t, inc, distance, order)

fig1 = plt.figure()
plt.plot(hp.sample_times,hp)
plt.plot(hc.sample_times,hc)
plt.axhline(y=0,color='k', linestyle='--')
plt.ticklabel_format(axis='y',style='sci',scilimits=(0,0),useMathText=True)
plt.legend(['hp','hx'])
plt.xlim([-0.01*M,0.01*M])
plt.savefig('refhphx_M_'+str(M)+'_q_'+str(q)+'.pdf')

psd = aLIGOZeroDetHighPower(len(hp) // 2 + 1, 1.0 / hp.duration, f_lower)

arrlen=50
eccs = np.linspace(1.07,1.5,arrlen)
impacts = np.linspace(40,200,arrlen)
eccs,impacts=np.meshgrid(eccs,impacts)

arr=np.zeros((arrlen,arrlen))
for i in tqdm(range(arrlen)):
    for j in range(arrlen):
        hp2,hc2=get_hyp_waveform(M,q,eccs[i][j], impacts[i][j],  delta_t, inc, distance, order)
        hp2 = hp2[:len(hp)] if len(hp) < len(hp2) else hp2
        hp2.resize(len(hp))

        m, idx = match(hp, hp2, psd=psd, low_frequency_cutoff=f_lower)
        arr[i][j]=m


import pickle

bigarr_file = open('file_M_'+str(M)+'_q_'+str(q)+'.pkl', 'ab')

pickle.dump(arr,bigarr_file)

bigarr_file.close()