
import scipy.stats as stats
import numpy as np
from astropy.modeling.models import Sersic1D
import matplotlib.pyplot as plt
print('Run')
mu = 0
variance = 1
sigma = np.sqrt(variance)
x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)

plt.figure()
y = stats.norm.pdf(x, mu, sigma)
y_fft = np.fft.fftshift(np.abs(np.fft.fft(y))) / np.sqrt(len(y))
plt.plot(x,y)
plt.plot(x,y_fft)
plt.savefig('./Gaucian.png', dpi=300, bbox_inches='tight')
