import numpy as np
from astropy.modeling.models import Sersic1D
from scipy.special import gammaincinv
import matplotlib.pyplot as plt

def sersic(r, r0, n, A0):

    I = A0 * np.exp(-gammaincinv(2 * n, 0.5) * ((r / r0) ** (1 / n) - 1) )

    return I

def I_R(Ie,R,Re,n):
    bn = 1.9992*n - 0.3271 #change to longer approx / or use b for n=4
    I = Ie*np.exp(-bn*((R/Re)**(1/n)-1))
    return I
