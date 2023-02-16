import numpy as np
from astropy.modeling.models import Sersic1D
from scipy.special import gammaincinv
import matplotlib.pyplot as plt

plt.close('all')

def sersic(r, r0, n, A0):

    # bn = 1.9992*n - 0.3271
    # bn = 2*n - 1./3. + (4. / 405.*n) + (46./2551.*n*n) + (131. / 1148175*n*n*n) - (2194697. / 30690717750*n*n*n*n)
    # I = A0*np.exp(-bn*((r/r0)**(1/n)-1))

    I = A0 * np.exp(-gammaincinv(2 * n, 0.5) * ((r / r0) ** (1 / n) - 1) )

    return I

r = np.linspace(0, 5, 1024)
r0 = 1
A0 = 1

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