import numpy as np
from astropy.modeling.models import Sersic1D
import matplotlib.pyplot as plt
print("Hello World")
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

#Sersic profile for different sersic profiles.
plt.savefig('./sersic_profile.png', dpi=300, bbox_inches='tight')

diff = (s1(r) - Ir[9])/s1(r)
mean_diff = np.mean(diff)
plt.figure()
plt.plot(r,diff)
plt.axhline(mean_diff, color = 'black', label = 'Mean')
plt.legend(loc = 'best')
plt.ylabel('Percentage Difference')
plt.xlabel('Radius')
plt.title('Difference between Astropy model and my model for a sersic index of 9')

plt.savefig('./diff_sersic_profile.png', dpi=300, bbox_inches='tight')

#example for the Andromeda galaxy(https://iopscience.iop.org/article/10.1088/0004-637X/739/1/20#apj398957t4)

R_A = 3.167 #arcseconds
Re_A = 0.69 #Kpsc0.69
Ie_A = 26.11
n_A = 1.80

r_A = np.arange(0,9,0.01) #in Kpc
plt.figure()
Ir_M31 = I_R(Ie_A ,r_A ,R_A ,n_A)
plt.plot(r_A,Ir_M31)

plt.xlabel('radius')
plt.ylabel('surface brightness')
plt.title('Sersic Profile for M31')

plt.savefig('./sersic_profile_M31.png', dpi=300, bbox_inches='tight')

#example for a galaxy twise as large as Andromeda

R_A = 3.167*2 #arcseconds
Re_A = 0.69 #Kpsc0.69
Ie_A = 26.11
n_A = 1.80

r_A = np.arange(0,9,0.01) #in Kpc

Ir_M31_2 = I_R(Ie_A ,r_A ,R_A ,n_A)

plt.plot(r_A,Ir_M31_2, label = '2 times RA')
plt.legend(loc = 'best')

plt.xlabel('radius')
plt.ylabel('surface brightness')
plt.title('Sersic Profile for M31_double_R_A')
#plt.savefig('./sersic_profile_double.png', dpi=300, bbox_inches='tight')

#example for a galaxy four times as large as Andromeda

R_A = 3.167*4 #arcseconds
Re_A = 0.69 #Kpsc0.69
Ie_A = 26.11
n_A = 1.80

r_A = np.arange(0,9,0.01) #in Kpc
#plt.figure()
plt.plot(r_A,I_R(Ie_A ,r_A ,R_A ,n_A), label = '4 times RA')
plt.legend(loc = 'best')

plt.xlabel('radius')
plt.ylabel('surface brightness')
plt.title('Sersic Profile for M31_quad_R_A')

plt.savefig('./sersic_profile_M31_2X_4X_RA.png', dpi=300, bbox_inches='tight')

#plotting multiple sersic indexes with the same radius

R_A = 3.167 #arcseconds
Re_A = 0.69 #Kpsc0.69
Ie_A = 26.11
n_A = 1.80

r_A = np.arange(0,9,0.01) #in Kpc
plt.figure()
plt.plot(r_A,I_R(Ie_A ,r_A ,R_A ,n_A), label = 'Sersic index = 1.8')
plt.plot(r_A,I_R(Ie_A ,r_A ,R_A ,2), label = 'Sersic index = 2')
plt.plot(r_A,I_R(Ie_A ,r_A ,R_A ,5), label = 'Sersic index = 5')
plt.plot(r_A,I_R(Ie_A ,r_A ,R_A ,7), label = 'Sersic index = 7')
plt.plot(r_A,I_R(Ie_A ,r_A ,R_A ,9), label = 'Sersic index = 9')

plt.legend(loc = 'best')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('radius')
plt.ylabel('surface brightness')
plt.title('Sersic Profile for M31, with differnt sersic indexes')

plt.savefig('./sersic_profile_M31_n.png', dpi=300, bbox_inches='tight')

#PLotting the forier transform of M31
#Fourier Transform of sersic profile
n=len(r_A)                     #number of radius points
dr=r_A[1]-r_A[0]                 # time interval

M31_fr = np.fft.fftshift(np.fft.fftfreq(n, dr))
M31_ft = np.fft.fftshift(np.fft.fft(Ir_M31))

#graph for sersic profile of M31
plt.figure(figsize=(20,12))
plt.subplot(321)
plt.plot(r_A,Ir_M31,'.-')
plt.xlabel('time (secs)')
plt.title('Sersic profile for M31')

plt.subplot(325)
plt.plot(M31_fr,np.abs(M31_ft), '.-')
plt.xlabel('freq (Hz)')
plt.title('spectrum, abs');
plt.savefig('./M31 sersic fourier transform.png', dpi=300, bbox_inches='tight')

#Plot a fourier transform of a galaxy twise the size of M31_fr
M31_ft_2 = np.fft.fftshift(np.fft.fft(Ir_M31_2))

#graph for sersic profile of M31
plt.figure(figsize=(20,12))
plt.subplot(321)
plt.plot(r_A,Ir_M31_2,'.-')
plt.xlabel('time (secs)')
plt.title('Sersic profile for a galaxy twise the radius of M31')

plt.subplot(325)
plt.plot(M31_fr,np.abs(M31_ft_2), '.-')
plt.xlabel('freq (Hz)')
plt.title('spectrum, abs');
plt.savefig('./M31x2 sersic fourier transform.png', dpi=300, bbox_inches='tight')