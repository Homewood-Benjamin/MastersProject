
import scipy.stats as stats
import numpy as np
from astropy.modeling.models import Sersic1D
import matplotlib.pyplot as plt
print('Run')
mu = 0
variance = 0.5
sigma = np.sqrt(variance)
#x = np.linspace(mu - 3*sigma, mu + 3*sigma, 500)
x = np.linspace(-3,3,500)

plt.figure()
y = stats.norm.pdf(x, mu, sigma)
y_fft = np.fft.fftshift(np.abs(np.fft.fft(y))) / np.sqrt(len(y))
#y_ifft = np.ifft.fftshift(np.abs(np.ifft.fft(y))) / np.sqrt(len(y))
plt.plot(x,y)
plt.plot(x,y_fft)
#plt.plot(x,y_ifft)
plt.savefig('./Gaucian.png', dpi=300, bbox_inches='tight')
