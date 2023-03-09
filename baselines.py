import numpy as np
import matplotlib.pyplot as plt

corrds = np.loadtxt('interferometer.csv', delimiter = ',', skiprows=1, usecols = (1,2))
long = corrds[:,0] * np.pi/180 #[rad]
lat = corrds[:,1] * np.pi/180 #[rad]

e_merlin = np.loadtxt('emerlin.csv', delimiter = ',', skiprows = 1)
emerlin_xx = e_merlin[:,0]
emerlin_yy = e_merlin[:,1]
emerlin_zz = e_merlin[:,2]

earth_rad = 6.3781*10**6 #[meters]

bin = np.array([2.00,2.25,2.50,2.75,3.00,3.25,3.50,3.75,4.00,4.25,4.50,4.75,5.00,5.25,5.50,5.75]) #bins for the hisogram

i = j = 0

#----------------------E_merlin-------------------------

def baseline (long,lat,height):
    long_diff = np.zeros((len(long),len(long)))

    for j in range(len(long)):
        for i in range(len(long)):
            long_diff[j,i] = long[j] - long[i]

    lat_diff = np.zeros((len(long),len(long)))

    for j in range(len(long)):
        for i in range(len(long)):
            lat_diff[j,i] = long[j] - long[i]

    height_diff = np.zeros((len(long),len(long)))

    for j in range(len(long)):
        for i in range(len(long)):
            height_diff[j,i] = long[j] - long[i]

    r = np.sqrt((long_diff)**2 + (lat_diff)**2 + (height_diff)**2)

    w_zeros = np.where(r!=0)

    r0 = np.log10(r[w_zeros[0],w_zeros[1]])

    return(r0)

e_merlin_baselines = baseline(emerlin_xx,emerlin_yy,emerlin_zz)

ang_emerlin = np.log10(206265*((1.4*10**9)/10**(e_merlin_baselines))) #[arcseconds]
#plt.savefig('./baselines_E-merlin2.png', dpi=300, bbox_inches='tight')

#----------------------------------VLA-----------------------------------
vla_posn = np.array(
    [
        [-1351.5783, 15411.6866, -36.9220],  # from jvla-m27_eall.ms
        [-125.0441, 1428.1956, -5.6652],
        [-8677.1357, -5860.2272, 19.2474],
        [6939.6252, -3239.1838, 42.7388],
        [-602.7972, 6880.1589, -17.0686],
        [-1654.8708, 18863.3822, -39.7787],
        [15538.4036, -7262.1508, 70.9228],
        [-401.3073, -270.5884, 2.2196],
        [4732.4674, -2208.1463, 28.9390],
        [-2642.5450, -1783.0279, 14.1914],
        [-410.9167, 4691.5365, -13.2065],
        [-1317.9043, -889.0530, 5.0181],
        [438.6360, -204.4722, 0.4958],
        [2888.9829, -1347.5757, 17.4316],
        [-14209.9382, -9604.7272, 20.8076],
        [-259.5562, 2963.4401, -8.5215],
        [-17388.6486, -11758.7053, 18.6633],
        [-6347.1614, -4285.2242, 11.9298],
        [-824.4586, 9407.5130, -22.1973],
        [12358.4512, -5773.1839, 46.3523],
        [1440.8421, -671.9128, 2.9230],
        [-4328.6309, -2921.5726, 15.4884],
        [-11302.3862, -7636.1216, 19.0827],
        [9487.4079, -4430.0736, 41.6457],
        [-38.1019, 434.7553, -1.3541],
        [-1074.4189, 12255.8562, -30.5254],
    ]
)

vla_xx = vla_posn[:, 0]
vla_yy = vla_posn[:, 1]
vla_zz = vla_posn[:, 2]

vla_baselines = baseline(vla_xx,vla_yy,vla_zz)
ang_size = np.log10(206265*((1.4*10**9)/10**(bin))) #[arcseconds]

#Define figure
fig, ax1 = plt.subplots()
bin = np.array([2.00,2.25,2.50,2.75,3.00,3.25,3.50,3.75,4.00,4.25,4.50,4.75,5.00,5.25,5.50,5.75])
#vla_freq = plt.hist(vla_baselines)
#vla_freq = vla_freq[0]

ax1.hist(vla_baselines, bins = bin)
ax1.hist(e_merlin_baselines, bins = bin)
plt.xlabel('Log10(baseline length [m])')
ax2 = ax1.twiny()
ax2.plot(ang_size,bin,'.', alpha = 0)
plt.xlabel('Log10 angular size at 1.4GHz')
plt.ylabel('Number of baselines')
plt.title('Baslines of E-merlin and VLA')

plt.savefig('./baselines_E-merlin_vla.png', dpi=300, bbox_inches='tight')

#-------------------for my data------------------------

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
