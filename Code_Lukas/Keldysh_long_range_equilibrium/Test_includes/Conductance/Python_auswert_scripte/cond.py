from numpy import *
import scipy.io as sio
import scipy.linalg
from matplotlib import rc
import matplotlib.pyplot as plt
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
import matplotlib.gridspec as gridspec

s="/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Conductance_zero_mag/conductance_zero_mag_N15_Vg0.250000_mu-1.475000_T0.010000_U00.550000.mat"	

str_int =empty([1500,60])
A = sio.loadmat(s)
y= A['y']
wf = A['wf'][0,:]
wbP = A['wbP'][0,:]
wbX = A['wbX'][0,:]
N = A['N'][0,0]

s=empty([len(wf)], dtype=complex)
p=empty([len(wbP)], dtype=complex)
x=empty([len(wbX)], dtype=complex)
d=empty([len(wbX)], dtype=complex)
for ind in range(len(wbP)):
 	p[ind] = y[0,5]['m'][0,ind]['m'][N,N]
for ind in range(len(wbX)):
 	x[ind] = y[0,7]['m'][0,ind]['m'][N,N]
 	d[ind] = y[0,10]['m'][0,ind]['m'][N,N]

#plt.plot(wf,real(b),color='blue')
#plt.plot(wf,imag(b),color='red')
plt.plot(wbP,imag(p),color='blue',marker='o')
plt.plot(wbX,imag(x),color='green',marker='o')
plt.plot(wbX,imag(d),color='red',marker='o')
plt.xlim([-5, 5])
#
#
#savename='cond.pdf' 
#plt.savefig(savename, format = 'pdf')
plt.show();
