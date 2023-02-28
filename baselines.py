import numpy as np
import matplotlib.pyplot as plt

corrds = np.loadtxt('interferometer.csv', delimiter = ',', skiprows=1, usecols = (1,2))
long = corrds[:,0] * np.pi/180 #[rad]
lat = corrds[:,1] * np.pi/180 #[rad]

earth_rad = 6.3781*10**6 #[meters]

i = j = 0

long_diff = np.zeros((len(long),len(long)))

for i in range(len(long)):
    long_diff[0,i] = long[0] - long[i]

for i in range(len(long)):
    long_diff[1,i] = long[1] - long[i]

for i in range(len(long)):
    long_diff[2,i] = long[2] - long[i]

for i in range(len(long)):
    long_diff[3,i] = long[3] - long[i]

for i in range(len(long)):
    long_diff[4,i] = long[4] - long[i]

for i in range(len(long)):
    long_diff[5,i] = long[5] - long[i]

for i in range(len(long)):
    long_diff[6,i] = long[6] - long[i]

for i in range(len(long)):
    long_diff[7,i] = long[7] - long[i]

long_dist_diff = 2*earth_rad*np.sin(long_diff/2)

lat_diff = np.zeros((len(long),len(long)))

for i in range(len(long)):
    lat_diff[0,i] = lat[0] - lat[i]

for i in range(len(long)):
    lat_diff[1,i] = lat[1] - lat[i]

for i in range(len(long)):
    lat_diff[2,i] = lat[2] - lat[i]

for i in range(len(long)):
    lat_diff[3,i] = lat[3] - lat[i]

for i in range(len(long)):
    lat_diff[4,i] = lat[4] - lat[i]

for i in range(len(long)):
    lat_diff[5,i] = lat[5] - lat[i]

for i in range(len(long)):
    lat_diff[6,i] = lat[6] - lat[i]

for i in range(len(long)):
    lat_diff[7,i] = lat[7] - lat[i]

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
