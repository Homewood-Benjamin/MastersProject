import numpy as np
import pylab as plt

print("Hello World")

import numpy as np
from astropy.modeling.models import Sersic1D
import matplotlib.pyplot as plt

#This is the inbuilt function based from astropy

plt.figure()
s1 = Sersic1D(amplitude=1, r_eff=1)
r=np.arange(0, 100, .01)

for n in range(1, 10):
     s1.n = n
     plt.loglog(r, s1(r), color=str(float(n) / 15))

plt.axis([1e-1, 30, 1e-2, 1e3])
plt.xlabel('log Radius')
plt.ylabel('log Surface Brightness')
plt.text(.25, 1.5, 'n=1')
plt.text(.25, 300, 'n=10')
plt.xticks([])
plt.yticks([])
plt.show()

plt.legend(loc = 'best')
plt.savefig('./sersic_profile_function.png', dpi=300, bbox_inches='tight')

#own made function

def I_R(Ie,R,Re,n):
  bn = 1.9992*n - 0.3271
  I = Ie*np.exp(-bn*((R/Re)**(1/n)-1))
  return I

Ir = np.zeros((10,len(r)))

for i in range(1,10):
  Ir[i,:]= (I_R(1,r,1,i))
  plt.plot(r,I_R(1,r,1,i))
  plt.xscale('log')
  plt.yscale('log')

plt.xlabel('log radius')
plt.ylabel('log surface brightness')
plt.axis([1e-1, 30, 1e-2, 1e3])

plt.savefig('./sersic_profile.png', dpi=300, bbox_inches='tight')

diff = (s1(r) - Ir[9])
plt.figure()
plt.semilogx(diff,r)
plt.ylim([0,1])
plt.ylabel('Difference')
plt.xlabel('Radius')
plt.title('Difference between inbuilt model and hand made model for a sersic index of 9')

plt.savefig('./diff_sersic_profile.png', dpi=300, bbox_inches='tight')

#example for the Andromeda galaxy(https://iopscience.iop.org/article/10.1088/0004-637X/739/1/20#apj398957t4)

R_A = 3.167 #arcseconds
Re_A = 0.69 #Kpsc0.69
Ie_A = 18.6
n_A = 1.80

r_A = np.arange(0,9,0.01) #in Kpc
plt.figure()
plt.plot(r_A,I_R(Ie_A,r_A,Re_A,n_A))
plt.xscale('log')
plt.yscale('log')

#ax = plt.gca()
#ax.invert_yaxis()

plt.xlabel('log radius Kpc')
plt.ylabel('log surface brightness')

plt.savefig('./sersic_profile_M31.png', dpi=300, bbox_inches='tight')


