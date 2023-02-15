
import scipy.stats as stats
import numpy as np
from astropy.modeling.models import Sersic1D
import matplotlib.pyplot as plt
print('Run')
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
plt.title('signal in time, sigma = 1.0 $\sigma$')

plt.subplot(325)
plt.plot(fr,np.abs(ft), '.-')
plt.xlabel('freq (Hz)')
plt.title('spectrum, abs');
plt.savefig('./Gaucian_ft.png', dpi=300, bbox_inches='tight')

#graph for 2 sigma
plt.figure(figsize=(20,12))
plt.subplot(321)
plt.plot(t,y2,'.-')
plt.xlabel('time (secs)')
plt.title('signal in time, for sigma = 0.1 $\sigma$')

plt.subplot(325)
plt.plot(fr,np.abs(ft2), '.-')
plt.xlabel('freq (Hz)')
plt.title('spectrum, abs');
plt.savefig('./Gaucian_ft2.png', dpi=300, bbox_inches='tight')
