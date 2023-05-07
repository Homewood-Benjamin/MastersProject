import numpy as np
from astropy.modeling.models import Sersic1D
import matplotlib.pyplot as plt
from scipy.special import gammaincinv
from sersic_func import sersic
from scipy import stats
print("Hello World")
#This is the inbuilt function based from astropy

plt.figure()
s1 = Sersic1D(amplitude=1, r_eff=1)
r=np.arange(0, 100, .01)

for n in range(1, 10):
     s1.n = n
     plt.loglog(r, s1(r), color=str(float(n) / 15))

Astro_s4 = Sersic1D(1,1,4)

plt.axis([1e-1, 30, 1e-2, 1e3])
plt.xlabel('log Radius')
plt.ylabel('log Surface Brightness')
plt.text(.25, 1.5, 'n=1')
plt.text(.25, 300, 'n=10')
plt.xticks([])
plt.yticks([])
plt.show()
plt.title('Astropy: sersic profile, showing surface brightness with radius')
plt.legend(loc = 'best')
plt.savefig('./sersic_profile_function.png', dpi=300, bbox_inches='tight')

r = np.linspace(0, 5, 1024)
r0 = 1
A0 = 1
## potting the difference plot ##

plt.figure(1, figsize=(2*4.5, 3.75))
for n in [1, 2, 3, 4]:

    sersic_ih = sersic(r, r0, n, A0)
    sersic_ast = Sersic1D(amplitude=A0, r_eff=r0, n=n)(r)

    plt.subplot(121)
    # offset by 10*n in the y direction to
    plt.plot(r, sersic_ih + 10*n, 'x', c='C{}'.format(n))
    plt.plot(r, sersic_ast + 10*n, '-', c='C{}'.format(n))
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('$r$')
    plt.ylabel('$I(r)$')
    plt.subplot(122)
    # also add a label to sanity check the max difference
    plt.plot(r, (sersic_ast - sersic_ih)/sersic_ast, '-',
             label='Max diff: {:.2e}'.format(np.max(sersic_ast - sersic_ih)))
    plt.xlabel('$r$')
    plt.ylabel('$\Delta I(r) / I_{\\rm astropy}(r)$')

plt.legend()
plt.subplots_adjust(wspace=0.3)
plt.show()

###

def I_R(Ie,R,Re,n):
    bn = 1.9992*n - 0.3271 #change to longer approx / or use b for n=4
    I = Ie*np.exp(-bn*((R/Re)**(1/n)-1))
    return I

Ir = np.zeros((10,len(r)))
plt.figure()
for i in range(1,10):
    Ir[i,:]= (sersic(r,1,i,1))
    plt.plot(r,sersic(r,1,i,1))
    plt.xscale('log')
    plt.yscale('log')

plt.legend(["n=1","n=2","n=3","n=4","n=5","n=6","n=7","n=8","n=9","n=10"])
plt.xlabel('log radius')
plt.ylabel('log surface brightness')
plt.title('Sersic profile, showing surface brightness with radius')
plt.axis([1e-1, 30, 1e-2, 1e3])

#Sersic profile for different sersic profiles.
plt.savefig('./sersic_profile.png', dpi=300, bbox_inches='tight')

diff = (s1(r) - Ir[9])/s1(r)
#diff = (Astro_s4(r) - Ir[4])/Astro_s4(r)
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

R_A = 33.73 #arcseconds
Re_A = 1.0 #Kpsc0.69
Ie_A = 26.11
n_A = 2.2

r_A = np.linspace(0,R_A,900) #in Kpc
plt.figure()
Ir_M31 = sersic(r_A,Re_A,n_A,Ie_A)

plt.plot(r_A,Ir_M31,label = 'r = RA')
#plt.gca().invert_yaxis()
plt.yscale('log')
plt.xscale('log')
plt.xlabel('radius [Kpsc]')
plt.ylabel('surface brightness')
plt.title('Sersic Profile for M31')

plt.savefig('./sersic_profile_M31.png', dpi=300, bbox_inches='tight')

#example for a galaxy twise as large as Andromeda

R_A = 33.73 #arcseconds
Re_A = 1.0*2 #Kpsc0.69
Ie_A = 26.11
n_A = 2.2

r_A = np.linspace(0,R_A,900) #in Kpc

Ir_M31_2 = sersic(r_A ,Re_A ,n_A,Ie_A)

plt.plot(r_A,Ir_M31_2, label = 'r = 2RA')
plt.legend(loc = 'best')

plt.xlabel('radius')
plt.ylabel('surface brightness')
plt.title('Sersic Profile for M31_double_R_A')
plt.savefig('./sersic_profile_double.png', dpi=300, bbox_inches='tight')

#example for a galaxy four times as large as Andromeda

R_A = 33.73 #arcseconds
Re_A = 1.0*4 #Kpsc0.69
Ie_A = 26.11
n_A = 2.2

r_A = np.linspace(0,R_A,900) #in Kpc
#plt.figure()
plt.plot(r_A,sersic(r_A ,Re_A ,n_A,Ie_A), label = 'r = 4RA')
plt.legend(loc = 'best')
#plt.xscale('log')
plt.yscale('log')
plt.xlabel('radius [Kpsc]')
plt.ylabel('surface brightness')
plt.title('Sersic Profile $R = Re,2Re,4Re$')

plt.savefig('./sersic_profile_M31_2X_4X_RA.png', dpi=300, bbox_inches='tight')

#plotting multiple sersic indexes with the same radius

R_A = 33.73 #arcseconds
Re_A = 1.0 #Kpsc0.69
Ie_A = 26.11
n_A = 2.2

Ir_n18 = sersic(r_A ,Re_A ,n_A,Ie_A )
Ir_n5 = sersic(r_A ,Re_A ,5,Ie_A)

r_A = np.linspace(0,R_A,900) #in Kpc

plt.figure()
plt.plot(r_A,sersic(r_A ,Re_A ,n_A,Ie_A), label = 'Sersic index = 2.2')
plt.plot(r_A,sersic(r_A ,Re_A ,3,Ie_A), label = 'Sersic index = 3')
plt.plot(r_A,sersic(r_A ,Re_A ,5,Ie_A), label = 'Sersic index = 5')
plt.plot(r_A,sersic(r_A ,Re_A ,7,Ie_A), label = 'Sersic index = 7')
plt.plot(r_A,sersic(r_A ,Re_A ,9,Ie_A), label = 'Sersic index = 9')

plt.legend(loc = 'best')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('radius')
plt.ylabel('surface brightness')
plt.title('Sersic Profile for M31, with differnt sersic indexes')

plt.savefig('./sersic_profile_M31_n.png', dpi=300, bbox_inches='tight')

##-------------------plotting a fourier transform of a gaussian----------------##
mu = 3
variance = 0.5
sigma = np.sqrt(variance)
x = np.arange(0,6,10)
y = stats.norm.pdf(x, mu, sigma)

fourierTransform = np.fft.fft(y)/(len(y))
fourierTransform = fourierTransform[range(int(len(y)/2))] # Exclude sampling frequency


n=1001                       #number of points
t = np.linspace(-3, 3, n )   #time
dt=t[1]-t[0]                 # time interval
print(f'sampling every  {dt:.3f} sec , so at {1/dt:.1f} Hz')

mu = 0
variance = 0.5
sigma = np.sqrt(variance)
y = stats.norm.pdf(t, mu, sigma)
y2 = stats.norm.pdf(t, mu, sigma*0.1)    # signal in time

fr= np.fft.fftshift(np.fft.fftfreq(n, dt))  # shift helps with sorting the frequencies for better plotting
ft=np.fft.fftshift(np.fft.fft(y))           # fftshift only necessary for plotting in sequence

ft2 = np.fft.fftshift(np.fft.fft(y2))

#graph for 1 sigma
plt.figure(figsize=(20,12))
plt.title('title')
plt.subplot(321)
plt.plot(t,y,'.-')
plt.xlabel('time (secs)')
plt.ylabel('y-axis')
plt.title('Signal in time, $\sigma$ = 1.0')

plt.subplot(325)
plt.plot(fr,np.abs(ft), '.-')
plt.xlabel('freq (Hz)')
plt.xlim(-25,25)
plt.ylabel('y-axis')
plt.title('Fourier Transform of Gaussian $\sigma$ = 1.0');
plt.savefig('./Gaucian_ft.png', dpi=300, bbox_inches='tight')

#graph for 2 sigma
plt.figure(figsize=(20,12))
plt.subplot(321)
plt.plot(t,y2,'.-')
plt.xlabel('time (secs)')
plt.ylabel('y-axis')
plt.title('Signal in time, for $\sigma$ = 0.1 $\sigma$')

plt.subplot(325)
plt.plot(fr,np.abs(ft2), '.-')
plt.xlabel('freq (Hz)')
plt.ylabel('y-axis')
plt.title('Fourier Transform of Gaussian $\sigma$ = 0.1');
plt.savefig('./Gaucian_ft2.png', dpi=300, bbox_inches='tight')

#---------------------------PLotting the forier transform of M31---------------------##
#Fourier Transform of sersic profile
n=len(r_A)                     #number of radius points
dr=r_A[1]-r_A[0]                 # time interval

M31_fr = np.fft.fftshift(np.fft.fftfreq(n, dr))
M31_ft = np.fft.fftshift(np.fft.fft(Ir_M31))

#graph for sersic profile of M31
plt.figure(figsize=(20,12))
plt.subplot(321)
plt.plot(r_A,Ir_M31,'.-')
plt.xlabel('Distance [Kpc]')
plt.title('Sersic profile for M31')

pos = np.argwhere(M31_fr>0)
M31_fr_pos = M31_fr[pos]
M31_ft_pos = M31_ft[pos]

plt.subplot(325)
plt.plot(M31_fr_pos,np.abs(M31_ft_pos), '.-')
plt.xlabel('$arcseconds^-1$')
plt.title('Sersic Profile $(Re = 1Re)$ Fourier Transform, abs');
plt.savefig('./M31_sersic_fourier_transform.png', dpi=300, bbox_inches='tight')

#Plot a fourier transform of a galaxy twise the size of M31_fr
M31_ft_2 = np.fft.fftshift(np.fft.fft(Ir_M31_2))

#graph for sersic profile of M31
plt.figure(figsize=(20,12))
plt.subplot(321)
plt.plot(r_A,Ir_M31_2,'.-')
plt.xlabel('Distance [Kpc]')
plt.title('Sersic profile for a galaxy twise the Effective Radius of M31')

pos_2re = np.argwhere(M31_fr>0)
M31_fr_pos = M31_fr[pos]
M31_ft_pos_2 = M31_ft_2[pos]

plt.subplot(325)
plt.plot(M31_fr_pos,np.abs(M31_ft_pos_2), '.-')
plt.xlabel('$arcseconds^-1$')
plt.title('Sersic Profile $(Re = 2Re)$ Fourier Transform, abs');

plt.savefig('./M31x2_sersic_fourier_transform.png', dpi=300, bbox_inches='tight')

#graph for sersic profile of M31 n=1.8
M31_ft_3 = np.fft.fftshift(np.fft.fft(Ir_n18))

plt.figure(figsize=(20,12))
plt.subplot(321)
plt.plot(r_A,Ir_n18,'.-')
plt.xlabel('Distance [Kpc]')
plt.title('Sersic profile for a galaxy M31 with n = 2')

pos_2n = np.argwhere(M31_fr>0)
M31_fr_pos = M31_fr[pos]
M31_ft_pos_3 = M31_ft_3[pos]

plt.subplot(325)
plt.plot(M31_fr_pos,np.abs(M31_ft_pos_3), '.-')
plt.xlabel('$arcseconds^-1$')
plt.title('Sersic Profile $(n = 2)$ Fourier Transform, abs');
plt.savefig('./M31n18_sersic_fourier_transform.png', dpi=300, bbox_inches='tight')

#graph for sersic profile of M31 n=5
M31_ft_4 = np.fft.fftshift(np.fft.fft(Ir_n5))

plt.figure(figsize=(20,12))
plt.subplot(321)
plt.plot(r_A,Ir_n5,'.-')
plt.xlabel('Distance [Kpc]')
plt.title('Sersic profile for a galaxy M31 with n = 5')

pos_5re = np.argwhere(M31_fr>0)
M31_fr_pos = M31_fr[pos]
M31_ft_pos_4 = M31_ft_4[pos]

plt.subplot(325)
plt.plot(M31_fr_pos,np.abs(M31_ft_pos_4), '.-')
plt.xlabel('$arcseconds^(-1)$')
plt.title('Sersic Profile $(n = 5)$ Fourier Transform, abs');
plt.savefig('./M31n5_sersic_fourier_transform.png', dpi=300, bbox_inches='tight')

#coppied from Git

import numpy as np
from astropy.modeling.models import Sersic1D
from scipy.special import gammaincinv
import matplotlib.pyplot as plt

plt.close('all')

r = np.linspace(0, 5, 1024)
r0 = 1
A0 = 1

plt.figure(1, figsize=(2*4.5, 3.75))
for n in [1, 2, 3, 4]:

    sersic_ih = sersic(A0, r, r0, n)
    sersic_ast = Sersic1D(amplitude=A0, r_eff=r0, n=n)(r)

    plt.subplot(121)
    # offset by 10*n in the y direction to
    plt.plot(r, sersic_ih + 10*n, 'x', c='C{}'.format(n))
    plt.plot(r, sersic_ast + 10*n, '-', c='C{}'.format(n))
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('$r$')
    plt.ylabel('$I(r)$')
    plt.subplot(122)
    # also add a label to sanity check the max difference
    plt.plot(r, (sersic_ast - sersic_ih)/sersic_ast, '-',
             label='Max diff: {:.2e}'.format(np.max(sersic_ast - sersic_ih)))
    plt.xlabel('$r$')
    plt.ylabel('$\Delta I(r) / I_{\\rm astropy}(r)$')

plt.legend()
plt.subplots_adjust(wspace=0.3)
plt.show()

##---------------------------discution plots----------------------##

plt.figure()
plt.subplot(231)
plt.title('Sersic Profile $(n = 2,5)$ Fourier Transform, abs')
plt.plot(M31_fr_pos,np.abs(M31_ft_pos_3), '.-', label = '$n=2$')
plt.xlabel('$arcseconds^(-1)$')
plt.legend(loc = 'best')
plt.yscale('log')
#plt.xscale('log')

plt.subplot(232)
plt.plot(M31_fr_pos,np.abs(M31_ft_pos_4), '.-', label = '$n=5$')
plt.xlabel('$arcseconds^(-1)$')
plt.legend(loc = 'best')
plt.yscale('log')
#plt.xscale('log')

plt.savefig('./M31n5_2_sersic_fourier_transform.png', dpi=300, bbox_inches='tight')

plt.figure()
plt.subplot(231)
plt.title('Sersic Profile $(n = 2,5)$')
plt.plot(r_A,Ir_n18, '.-', label = '$n=2$')
plt.xlabel('Distance [Kpsc]')
plt.legend(loc = 'best')
plt.yscale('log')
#plt.xscale('log')

plt.subplot(232)
plt.plot(r_A,Ir_n5, '.-', label = '$n=5$')
plt.xlabel('Distance [Kpsc]')
plt.legend(loc = 'best')
plt.yscale('log')
#plt.xscale('log')

plt.savefig('./M31n5_2_sersic.png', dpi=300, bbox_inches='tight')

plt.figure()
plt.plot(M31_fr_pos,np.abs(M31_ft_pos), '.-', label = '$Re = 1Re$')
plt.xlabel('$arcseconds^-1$')
plt.legend(loc = 'best')
plt.savefig('./M31_1RE_sersic_fourier_transform.png', dpi=300, bbox_inches='tight')

plt.plot(M31_fr_pos,np.abs(M31_ft_pos_2), '.-', label = '$Re = 2Re$')
plt.xlabel('$arcseconds^-1$')
plt.legend(loc = 'best')
plt.title('Sersic Profile $(Re = Re & 2Re)$ Fourier Transform, abs');

plt.savefig('./M31_2RE_sersic_fourier_transform.png', dpi=300, bbox_inches='tight')
