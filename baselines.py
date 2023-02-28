import numpy as np
import matplotlib.pyplot as plt

corrds = np.loadtxt('interferometer.csv', delimiter = ',', skiprows=1, usecols = (1,2))
long = corrds[:,0] * np.pi/180 #[rad]
lat = corrds[:,1] * np.pi/180 #[rad]

earth_rad = 6.3781*10**6 #[meters]

i = j = 0

long_diff = np.zeros((len(long),len(long)))

for j in range(len(long)):
    for i in range(len(long)):
        long_diff[j,i] = long[j] - long[i]

long_dist_diff = 2*earth_rad*np.sin(long_diff/2)

lat_diff = np.zeros((len(long),len(long)))

for j in range(len(long)):
    for i in range(len(long)):
        lat_diff[j,i] = long[j] - long[i]

lat_dist_diff = 2*earth_rad*np.sin(lat_diff/2)

r = np.sqrt((long_dist_diff)**2 + (lat_dist_diff)**2)

#plot baselines as a histogram

w_zeros = np.where(r!=0)

r0 = np.log10(r[w_zeros[0],w_zeros[1]])

plt.figure()
plt.hist(r0, bins = [2.00,2.25,2.50,2.75,3.00,3.25,3.50,3.75,4.00,4.25,4.50,4.75,5.00,5.25,5.50,5.75])
plt.xlabel('Log10(baseline length [m])')
plt.ylabel('Number of baselines)')
plt.title('Baslines of E-merlin')

plt.savefig('./baselines_E-merlin.png', dpi=300, bbox_inches='tight')
