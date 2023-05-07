import numpy as np
import matplotlib.pyplot as plt
from sersic_func import sersic

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
lam = (3.0*10**8)/(1.4*10**9) #[Hz]
#Define figure
fig, ax1 = plt.subplots()

bin = np.array([2.00,2.25,2.50,2.75,3.00,3.25,3.50,3.75,4.00,4.25,4.50,4.75,5.00,5.25,5.50,5.75])

from astropy import units

ang_size = ((lam / 10**bin)*units.rad).to(units.arcsec)
#vla_freq = plt.hist(vla_baselines)
#vla_freq = vla_freq[0]

R_A = 3.167 #arcseconds
Re_A = 0.69 #Kpsc0.69
Ie_A = 26.11
n_A = 1.80

r_A = np.arange(0,9,0.01) #in Kpc

Ir_M31 = sersic(r_A,R_A,n_A,Ie_A)

n=len(r_A)                     #number of radius points
dr=r_A[1]-r_A[0]                 # time interval

M31_fr = np.fft.fftshift(np.fft.fftfreq(n, dr))
M31_ft = np.fft.fftshift(np.fft.fft(Ir_M31))

#graph for sersic profile of M31
#plt.plot(M31_fr,np.abs(M31_ft), '.-')

ax1.hist(vla_baselines, bins = bin, label = 'VLA')
ax1.hist(e_merlin_baselines, bins = bin, label = 'E-merlin', color = "tab:orange")
plt.legend(loc = 'best')
ax2 = ax1.twiny()
ax2.plot(ang_size,bin,'.', alpha = 0)
plt.xscale('log')
ax2.invert_xaxis()
ax2.set_xlabel('Log10 angular size at 1.4GHz')
ax1.set_xlabel('Log10(baseline Length [m])')
ax1.set_ylabel('Number of baselines')
plt.title('Baslines of E-merlin and VLA')

plt.savefig('./baselines_E-merlin_vla.png', dpi=300, bbox_inches='tight')
#
#------------------combining the baselines plots with sersic profiles----------------
##INFO FOR SERSIC PROFILE/Fourier Transform##
R_A = 33.37 #3.167 #arcseconds
Re_A = 1.0 # 0.69 #Kpsc0.69
Ie_A = 26.11
n_A = 2 #1.80

r_A = np.linspace(0,R_A,900) #in Kpc

Ir_M31 = sersic(r_A,R_A,n_A,Ie_A)
Ir_M31_2 = sersic(r_A ,Re_A*2 ,n_A,Ie_A)
Ir_M31_3 = sersic(r_A ,Re_A*3 ,n_A,Ie_A)
Ir_M31_4 = sersic(r_A ,Re_A*4 ,n_A,Ie_A)
Ir_M31_5 = sersic(r_A ,Re_A*5 ,n_A,Ie_A)

n=len(r_A)                     #number of radius points
dr=r_A[1]-r_A[0]                 # time interval

M31_fr = np.fft.fftshift(np.fft.fftfreq(n, dr))
M31_ft = np.fft.fftshift(np.fft.fft(Ir_M31))
M31_ft_2 = np.fft.fftshift(np.fft.fft(Ir_M31_2))
M31_ft_3 = np.fft.fftshift(np.fft.fft(Ir_M31_3))
M31_ft_4 = np.fft.fftshift(np.fft.fft(Ir_M31_4))
M31_ft_5 = np.fft.fftshift(np.fft.fft(Ir_M31_5))

pos = np.argwhere(M31_fr>0)
M31_fr_pos = M31_fr[pos]
M31_ft_pos = M31_ft[pos]/max(M31_ft[pos]) #Normalised

pos_2 = np.argwhere(M31_fr>0)
M31_fr_pos_2 = M31_fr[pos_2]
M31_ft_pos_2 = M31_ft_2[pos]/max(M31_ft_2[pos_2]) #Normalised

pos_3 = np.argwhere(M31_fr>0)
M31_fr_pos_3 = M31_fr[pos_3]
M31_ft_pos_3 = M31_ft_3[pos_3]/max(M31_ft_3[pos_3]) #Normalised

pos_4 = np.argwhere(M31_fr>0)
M31_fr_pos_4 = M31_fr[pos_4]
M31_ft_pos_4 = M31_ft_4[pos_4]/max(M31_ft_4[pos_4]) #Normalised

pos_5 = np.argwhere(M31_fr>0)
M31_fr_pos_5 = M31_fr[pos_5]
M31_ft_pos_5 = M31_ft_5[pos_5]/max(M31_ft_5[pos_5]) #Normalised

fig, ax1 = plt.subplots()

ang_size_vla_baselines = (lam / 10**vla_baselines)*((3600 * 180)/np.pi)
ang_size_e_merlin_baselines = (lam / 10**e_merlin_baselines)*((3600 * 180)/np.pi)
ang_size_bin = (lam / 10**bin)*((3600 * 180)/np.pi)
ax1.hist(ang_size_vla_baselines, bins = ang_size_bin[::-1], label = 'VLA', color = 'C0')
ax1.hist(ang_size_e_merlin_baselines,bins = ang_size_bin[::-1], label = 'E-merlin',color = 'C1')


plt.xscale('log')


ax1.plot(M31_fr_pos,175*np.abs(M31_ft_pos), label = 'r = RA' ,color = 'black')
ax1.plot(M31_fr_pos_2,175*np.abs(M31_ft_pos_2), label = 'r = 2RA' ,color = 'blue')
ax1.plot(M31_fr_pos_3,175*np.abs(M31_ft_pos_3), label = 'r = 3RA' ,color = 'green')
ax1.plot(M31_fr_pos_4,175*np.abs(M31_ft_pos_4), label = 'r = 4RA', color = 'red')
ax1.plot(M31_fr_pos_5,175*np.abs(M31_ft_pos_5), label = 'r = 5RA' ,color = 'orange')

#plt.xlim(0,10)
#plt.ylim(0,175)
plt.legend(loc = 'best')
plt.xlabel('Arcseconds')
plt.ylabel('Normalised Brightness')

plt.savefig('./fourier_baselines.png', dpi=300, bbox_inches='tight')

##-----------------changing the sersic index-------------

Ir_M31 = sersic(r_A,Re_A,1,Ie_A)
Ir_M31_n2 = sersic(r_A ,Re_A ,2,Ie_A)
Ir_M31_n3 = sersic(r_A ,Re_A ,3,Ie_A)
Ir_M31_n4 = sersic(r_A ,Re_A ,4,Ie_A)
Ir_M31_n5 = sersic(r_A ,Re_A ,5,Ie_A)

M31_fr = np.fft.fftshift(np.fft.fftfreq(n, dr))
M31_ft_n = np.fft.fftshift(np.fft.fft(Ir_M31))
M31_ft_2n = np.fft.fftshift(np.fft.fft(Ir_M31_n2))
M31_ft_3n = np.fft.fftshift(np.fft.fft(Ir_M31_n3))
M31_ft_4n = np.fft.fftshift(np.fft.fft(Ir_M31_n4))
M31_ft_5n = np.fft.fftshift(np.fft.fft(Ir_M31_n5))

pos_n = np.argwhere(M31_fr>0)
M31_fr_pos_n = M31_fr[pos]
M31_ft_pos_n = M31_ft_n[pos_n]/max(M31_ft_n[pos_n]) #Normalised

pos_2n = np.argwhere(M31_fr>0)
M31_fr_pos_2n = M31_fr[pos_2n]
M31_ft_pos_2n = M31_ft_2n[pos]/max(M31_ft_2n[pos_2n]) #Normalised

pos_3n = np.argwhere(M31_fr>0)
M31_fr_pos_3n = M31_fr[pos_3n]
M31_ft_pos_3n = M31_ft_3n[pos_3n]/max(M31_ft_3n[pos_3n]) #Normalised

pos_4n = np.argwhere(M31_fr>0)
M31_fr_pos_4n = M31_fr[pos_4n]
M31_ft_pos_4n = M31_ft_4n[pos_4n]/max(M31_ft_4n[pos_4n]) #Normalised

pos_5n = np.argwhere(M31_fr>0)
M31_fr_pos_5n = M31_fr[pos_5n]
M31_ft_pos_5n = M31_ft_5n[pos_5n]/max(M31_ft_5n[pos_5n]) #Normalised


plt.figure()
fig, ax1 = plt.subplots()
ax1.hist(ang_size_vla_baselines, bins = ang_size_bin[::-1], label = 'VLA', color = 'C0')
ax1.hist(ang_size_e_merlin_baselines,bins = ang_size_bin[::-1], label = 'E-merlin',color = 'C1')

plt.xscale('log')


ax1.plot(M31_fr_pos_n,175*np.abs(M31_ft_pos), label = 'n = 1.8' ,color = 'black')
ax1.plot(M31_fr_pos_2n,175*np.abs(M31_ft_pos_2n), label = 'n = 2' ,color = 'blue')
ax1.plot(M31_fr_pos_3n,175*np.abs(M31_ft_pos_3n), label = 'n = 3' ,color = 'green')
ax1.plot(M31_fr_pos_4n,175*np.abs(M31_ft_pos_4n), label = 'n = 4', color = 'red')
ax1.plot(M31_fr_pos_5n,175*np.abs(M31_ft_pos_5n), label = 'n = 5' ,color = 'orange')

#plt.xlim(0,50)
#plt.ylim(0,175)
plt.legend(loc = 'best')
plt.xlabel('Arcseconds')
plt.ylabel('Normalised Brightness')

plt.savefig('./fourier_baselines_n.png', dpi=300, bbox_inches='tight')

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
plt.ylabel('Number of baselines')
plt.title('Baslines of E-merlin')

plt.savefig('./baselines_E-merlin_mydata.png', dpi=300, bbox_inches='tight')

#----------------------------------------------------------------------------##
#in_bin = np.zeros((len(ang_size_bin)))

in_bin0 = np.argwhere((M31_fr_pos<ang_size_bin[0]) & (M31_fr_pos>ang_size_bin[1]))
in_bin1 = np.argwhere((M31_fr_pos<ang_size_bin[1]) & (M31_fr_pos>ang_size_bin[2]))
in_bin2 = np.argwhere((M31_fr_pos<ang_size_bin[2]) & (M31_fr_pos>ang_size_bin[3]))
in_bin3 = np.argwhere((M31_fr_pos<ang_size_bin[3]) & (M31_fr_pos>ang_size_bin[4]))
in_bin4 = np.argwhere((M31_fr_pos<ang_size_bin[4]) & (M31_fr_pos>ang_size_bin[5]))
in_bin5 = np.argwhere((M31_fr_pos<ang_size_bin[5]) & (M31_fr_pos>ang_size_bin[6]))
in_bin6 = np.argwhere((M31_fr_pos<ang_size_bin[6]) & (M31_fr_pos>ang_size_bin[7]))
in_bin7 = np.argwhere((M31_fr_pos<ang_size_bin[7]) & (M31_fr_pos>ang_size_bin[8]))
in_bin8 = np.argwhere((M31_fr_pos<ang_size_bin[8]) & (M31_fr_pos>ang_size_bin[9]))
in_bin9 = np.argwhere((M31_fr_pos<ang_size_bin[9]) & (M31_fr_pos>ang_size_bin[10]))
in_bin10 = np.argwhere((M31_fr_pos<ang_size_bin[10]) & (M31_fr_pos>ang_size_bin[11]))
in_bin11 = np.argwhere((M31_fr_pos<ang_size_bin[11]) & (M31_fr_pos>ang_size_bin[12]))
in_bin12 = np.argwhere((M31_fr_pos<ang_size_bin[12]) & (M31_fr_pos>ang_size_bin[13]))
in_bin13 = np.argwhere((M31_fr_pos<ang_size_bin[13]) & (M31_fr_pos>ang_size_bin[14]))
in_bin14 = np.argwhere((M31_fr_pos<ang_size_bin[14]) & (M31_fr_pos>ang_size_bin[15]))
in_bin15 = np.argwhere(M31_fr_pos<ang_size_bin[15])

inbin = np.array([in_bin0,in_bin1,in_bin2,in_bin3,in_bin4,in_bin5,in_bin6,in_bin7,in_bin8,in_bin9,in_bin10,in_bin11,in_bin12,in_bin13,in_bin14])

vla_height = ax1.hist(ang_size_vla_baselines, bins = ang_size_bin[::-1], label = 'VLA', color = 'C0')[0]
def error(N):
    aprox_error = 420/(0.805*np.sqrt(2*N*(N-1))) #add a time of 24hrs in sec
    return aprox_error
vla_error_aprox = np.zeros(len(vla_height))

for i in range(len(vla_height)):
    vla_error_aprox[i] = error(vla_height[i])



w_M31_fr_pos0 = vla_error_aprox[0]+M31_ft_pos[inbin[0][:,0]] #y-axis
w_M31_fr_pos1 = vla_error_aprox[1]+M31_ft_pos[inbin[1][:,0]] #y-axis
w_M31_fr_pos2 = vla_error_aprox[2]+M31_ft_pos[inbin[2][:,0]] #y-axis
w_M31_fr_pos3 = vla_error_aprox[3]+M31_ft_pos[inbin[3][:,0]] #y-axis
w_M31_fr_pos4 = vla_error_aprox[4]+M31_ft_pos[inbin[4][:,0]] #y-axis
w_M31_fr_pos5 = vla_error_aprox[5]+M31_ft_pos[inbin[5][:,0]] #y-axis
w_M31_fr_pos6 = vla_error_aprox[6]+M31_ft_pos[inbin[6][:,0]] #y-axis
w_M31_fr_pos7 = vla_error_aprox[7]+M31_ft_pos[inbin[7][:,0]] #y-axis
w_M31_fr_pos8 = vla_error_aprox[8]+M31_ft_pos[inbin[8][:,0]] #y-axis
w_M31_fr_pos9 = vla_error_aprox[9]+M31_ft_pos[inbin[9][:,0]] #y-axis
w_M31_fr_pos10 = vla_error_aprox[10]+M31_ft_pos[inbin[10][:,0]] #y-axis
w_M31_fr_pos11 = vla_error_aprox[11]+M31_ft_pos[inbin[11][:,0]] #y-axis
w_M31_fr_pos12 = vla_error_aprox[12]+M31_ft_pos[inbin[12][:,0]] #y-axis
w_M31_fr_pos13 = vla_error_aprox[13]+M31_ft_pos[inbin[13][:,0]] #y-axis
w_M31_fr_pos14 = vla_error_aprox[14]+M31_ft_pos[inbin[14][:,0]] #y-axis

w_M31_ft_pos = np.concatenate((w_M31_fr_pos0,w_M31_fr_pos1,w_M31_fr_pos2,w_M31_fr_pos3,w_M31_fr_pos4,w_M31_fr_pos5,w_M31_fr_pos6,w_M31_fr_pos7,w_M31_fr_pos8,w_M31_fr_pos9,w_M31_fr_pos10,w_M31_fr_pos11,w_M31_fr_pos12,w_M31_fr_pos13,w_M31_fr_pos14))

#w_M31_fr_pos = np.zeros(len(vla_height))
#for i in range(len(vla_height)):
    #w_M31_fr_pos[i] = vla_error_aprox[i]+M31_ft_pos[inbin[i][:,0]] #y-axis



plt.figure()
plt.xscale('log')
#plt.plot(M31_fr_pos,w_M31_ft_pos)
plt.xlabel('arcseconds')
plt.ylabel('Intensity')
plt.savefig('./weighted_baselines_VLA.png', dpi=300, bbox_inches='tight')

##---------attampt at 'weighting' the sersic profile with errors---------------------------------------------------##

#ang_size_bin = ang_size_bin[::-1]

mid_point = 10**(bin+(0.25/2))
mid_point_ang = (lam / mid_point)*((3600 * 180)/np.pi)
mid_point_ang = mid_point_ang[::-1]


bin_range = np.zeros((len(bin),2))

for i in range(len(bin)):
    bin_range[i,0] = bin[i]
for i in range(len(bin)-1):
    bin_range[i,1] = bin[i+1]

bin_range = 10**(bin_range)
bin_range[-1,-1] = 0

bin_range_ang = (lam / 10**bin_range)*((3600 * 180)/np.pi)

#M31_fr_pos_ang = (lam / 10**M31_fr_pos)*((3600 * 180)/np.pi)

sersic_X_b0 = np.argwhere((M31_fr_pos>ang_size_bin[::-1][0]) & (M31_fr_pos<ang_size_bin[::-1][1]))
sersic_X_b1 = np.argwhere((M31_fr_pos>ang_size_bin[::-1][1]) & (M31_fr_pos<ang_size_bin[::-1][2]))
sersic_X_b2 = np.argwhere((M31_fr_pos>ang_size_bin[::-1][2]) & (M31_fr_pos<ang_size_bin[::-1][3]))
sersic_X_b3 = np.argwhere((M31_fr_pos>ang_size_bin[::-1][3]) & (M31_fr_pos<ang_size_bin[::-1][4]))
sersic_X_b4 = np.argwhere((M31_fr_pos>ang_size_bin[::-1][4]) & (M31_fr_pos<ang_size_bin[::-1][5]))
sersic_X_b5 = np.argwhere((M31_fr_pos>ang_size_bin[::-1][5]) & (M31_fr_pos<ang_size_bin[::-1][6]))
sersic_X_b6 = np.argwhere((M31_fr_pos>ang_size_bin[::-1][6]) & (M31_fr_pos<ang_size_bin[::-1][7]))
sersic_X_b7 = np.argwhere((M31_fr_pos>ang_size_bin[::-1][7]) & (M31_fr_pos<ang_size_bin[::-1][8]))
sersic_X_b8 = np.argwhere((M31_fr_pos>ang_size_bin[::-1][8]) & (M31_fr_pos<ang_size_bin[::-1][9]))
sersic_X_b9 = np.argwhere((M31_fr_pos>ang_size_bin[::-1][9]) & (M31_fr_pos<ang_size_bin[::-1][10]))
sersic_X_b10 = np.argwhere((M31_fr_pos>ang_size_bin[::-1][10]) & (M31_fr_pos<ang_size_bin[::-1][11]))
sersic_X_b11 = np.argwhere((M31_fr_pos>ang_size_bin[::-1][11]) & (M31_fr_pos<ang_size_bin[::-1][12]))
sersic_X_b12 = np.argwhere((M31_fr_pos>ang_size_bin[::-1][12]) & (M31_fr_pos<ang_size_bin[::-1][13]))
sersic_X_b13 = np.argwhere((M31_fr_pos>ang_size_bin[::-1][13]) & (M31_fr_pos<ang_size_bin[::-1][14]))
sersic_X_b14 = np.argwhere((M31_fr_pos>ang_size_bin[::-1][14]) & (M31_fr_pos<ang_size_bin[::-1][15]))


in_bin_0 =  M31_fr_pos[sersic_X_b0][:,0]
in_bin_1 =  M31_fr_pos[sersic_X_b1][:,0]
in_bin_2 =  M31_fr_pos[sersic_X_b2][:,0]
in_bin_3 =  M31_fr_pos[sersic_X_b3][:,0]
in_bin_4 =  M31_fr_pos[sersic_X_b4][:,0]
in_bin_5 =  M31_fr_pos[sersic_X_b5][:,0]
in_bin_6 =  M31_fr_pos[sersic_X_b6][:,0]
in_bin_7 =  M31_fr_pos[sersic_X_b7][:,0]
in_bin_8 =  M31_fr_pos[sersic_X_b8][:,0]
in_bin_9 =  M31_fr_pos[sersic_X_b9][:,0]
in_bin_10 =  M31_fr_pos[sersic_X_b10][:,0]
in_bin_11 =  M31_fr_pos[sersic_X_b11][:,0]
in_bin_12 =  M31_fr_pos[sersic_X_b12][:,0]
in_bin_13 =  M31_fr_pos[sersic_X_b13][:,0]
in_bin_14 =  M31_fr_pos[sersic_X_b14][:,0]

I_in_bin_0 = np.abs(M31_ft_pos[sersic_X_b0][:,0])
I_in_bin_1 =  np.abs(M31_ft_pos[sersic_X_b1][:,0])
I_in_bin_2 =  np.abs(M31_ft_pos[sersic_X_b2][:,0])
I_in_bin_3 =  np.abs(M31_ft_pos[sersic_X_b3][:,0])
I_in_bin_4 =  np.abs(M31_ft_pos[sersic_X_b4][:,0])
I_in_bin_5 =  np.abs(M31_ft_pos[sersic_X_b5][:,0])
I_in_bin_6 =  np.abs(M31_ft_pos[sersic_X_b6][:,0])
I_in_bin_7 =  np.abs(M31_ft_pos[sersic_X_b7][:,0])
I_in_bin_8 =  np.abs(M31_ft_pos[sersic_X_b8][:,0])
I_in_bin_9 =  np.abs(M31_ft_pos[sersic_X_b9][:,0])
I_in_bin_10 =  np.abs(M31_ft_pos[sersic_X_b10][:,0])
I_in_bin_11 =  np.abs(M31_ft_pos[sersic_X_b11][:,0])
I_in_bin_12 =  np.abs(M31_ft_pos[sersic_X_b12][:,0])
I_in_bin_13 =  np.abs(M31_ft_pos[sersic_X_b13][:,0])
I_in_bin_14 =  np.abs(M31_ft_pos[sersic_X_b14][:,0])

I_in_bin = np.concatenate((I_in_bin_0,I_in_bin_1,I_in_bin_2,I_in_bin_3,I_in_bin_4,I_in_bin_5,I_in_bin_6,I_in_bin_7,I_in_bin_8,I_in_bin_9,I_in_bin_10,I_in_bin_11,I_in_bin_12,I_in_bin_13,I_in_bin_14))

index_b0 = len(I_in_bin_0)
index_b1 = index_b0 + len(I_in_bin_1)
index_b2 = index_b1 + len(I_in_bin_2)
index_b3 = index_b2 + len(I_in_bin_3)
index_b4 = index_b3 + len(I_in_bin_4)
index_b5 = index_b4 + len(I_in_bin_5)
index_b6 = index_b5 + len(I_in_bin_6)
index_b7 = index_b6 + len(I_in_bin_7)
index_b8 = index_b7 + len(I_in_bin_8)
index_b9 = index_b8 + len(I_in_bin_9)
index_b10 = index_b9 + len(I_in_bin_10)
index_b11 = index_b10 + len(I_in_bin_11)
index_b12 = index_b11 + len(I_in_bin_12)
index_b13 = index_b12 + len(I_in_bin_13)
index_b14 = index_b13 + len(I_in_bin_14)

def dela_I(SEFD,nc,n_p,N,t,v):
    dI = SEFD/(nc*np.sqrt(n_p*N*(N-1)*t*v))
    return dI

SEFD = 420 #Jy
nc = 0.93 #quantum reciver effiency
n_p = 2 #number of polarisations
N = vla_height
t = 86400 #arcseconds
v = 1 #GHz

VLA_dI = dela_I(SEFD,nc,n_p,N,t,v)
#VLA_dI[0] = VLA_dI[1] = VLA_dI[2] = 0

'''rand_err0 = np.random.normal(175*I_in_bin_0,(VLA_dI[0]/2))
rand_err1 = np.random.normal(175*I_in_bin_1,(VLA_dI[1]/2))
rand_err2 = np.random.normal(175*I_in_bin_2,(VLA_dI[2]/2))
rand_err3 = np.random.normal(175*I_in_bin_3,(VLA_dI[3]/2))
rand_err4 = np.random.normal(175*I_in_bin_4,(VLA_dI[4]/2))
rand_err5 = np.random.normal(175*I_in_bin_5,(VLA_dI[5]/2))
rand_err6 = np.random.normal(175*I_in_bin_6,(VLA_dI[6]/2))
rand_err7 = np.random.normal(175*I_in_bin_7,(VLA_dI[7]/2))
rand_err8 = np.random.normal(175*I_in_bin_8,(VLA_dI[8]/2))
rand_err9 = np.random.normal(175*I_in_bin_9,(VLA_dI[9]/2))
rand_err10 = np.random.normal(175*I_in_bin_10,(VLA_dI[10]/2))
rand_err11 = np.random.normal(175*I_in_bin_11,(VLA_dI[11]/2))
rand_err12 = np.random.normal(175*I_in_bin_12,(VLA_dI[12]/2))
rand_err13 = np.random.normal(175*I_in_bin_13,(VLA_dI[13]/2))
rand_err14 = np.random.normal(175*I_in_bin_14,(VLA_dI[14]/2))'''

#rand_err0 = np.random.normal(0,(VLA_dI[0]))
#rand_err1 = np.random.normal(0,(VLA_dI[1]))
#rand_err2 = np.random.normal(0,(VLA_dI[2]))
rand_err3 = np.random.normal(0,(VLA_dI[3]))
rand_err4 = np.random.normal(0,(VLA_dI[4]))
rand_err5 = np.random.normal(0,(VLA_dI[5]))
rand_err6 = np.random.normal(0,(VLA_dI[6]))
rand_err7 = np.random.normal(0,(VLA_dI[7]))
rand_err8 = np.random.normal(0,(VLA_dI[8]))
rand_err9 = np.random.normal(0,(VLA_dI[9]))
rand_err10 = np.random.normal(0,(VLA_dI[10]))
rand_err11 = np.random.normal(0,(VLA_dI[11]))
rand_err12 = np.random.normal(0,(VLA_dI[12]))
rand_err13 = np.random.normal(0,(VLA_dI[13]))
rand_err14 = np.random.normal(0,(VLA_dI[14]))

rand_err = np.array((rand_err3,rand_err4,rand_err5,rand_err6,rand_err7,rand_err8,rand_err9,rand_err10,rand_err11,rand_err12,rand_err13,rand_err14))
#rand_err = rand_err[:,0]

plt.figure()
plt.plot(np.abs(rand_err))
#plt.xscale('log')
plt.savefig('./rand_err.png', dpi=300, bbox_inches='tight')


#sim_data0 = I_in_bin_0 + rand_err0
#sim_data1 = I_in_bin_1 + rand_err1
#sim_data2 = I_in_bin_2 + rand_err2
sim_data3 = I_in_bin_3 + rand_err3
sim_data4 = I_in_bin_4 + rand_err4
sim_data5 = I_in_bin_5 + rand_err5
sim_data6 = I_in_bin_6 + rand_err6
sim_data7 = I_in_bin_7 + rand_err7
sim_data8 = I_in_bin_8 + rand_err8
sim_data9 = I_in_bin_9 + rand_err9
sim_data10 = I_in_bin_10 + rand_err10
sim_data11 = I_in_bin_11 + rand_err11
sim_data12 = I_in_bin_12 + rand_err12
sim_data13 = I_in_bin_13 + rand_err13
sim_data14 = I_in_bin_14 + rand_err14

rand_err_3 = np.array([rand_err3]*len(I_in_bin_3))
rand_err_4 = np.array([rand_err4]*len(I_in_bin_4))
rand_err_5 = np.array([rand_err5]*len(I_in_bin_5))
rand_err_6 = np.array([rand_err6]*len(I_in_bin_6))
rand_err_7 = np.array([rand_err7]*len(I_in_bin_7))
rand_err_8 = np.array([rand_err8]*len(I_in_bin_8))
rand_err_9 = np.array([rand_err9]*len(I_in_bin_9))
rand_err_10 = np.array([rand_err10]*len(I_in_bin_10))
rand_err_11= np.array([rand_err11]*len(I_in_bin_11))
rand_err_12 = np.array([rand_err12]*len(I_in_bin_12))
rand_err_13 = np.array([rand_err13]*len(I_in_bin_13))
rand_err_14 = np.array([rand_err14]*len(I_in_bin_14))

randerr = np.concatenate((rand_err_3,rand_err_4,rand_err_5,rand_err_6,rand_err_7,rand_err_8,rand_err_9,rand_err_10,rand_err_11,rand_err_12,rand_err_13,rand_err_14))

'''sim_data0 = I_in_bin_0 + VLA_dI[0]
sim_data1 = I_in_bin_1 + VLA_dI[1]
sim_data2 = I_in_bin_2 + VLA_dI[2]
sim_data3 = I_in_bin_3 + VLA_dI[3]
sim_data4 = I_in_bin_4 + VLA_dI[4]
sim_data5 = I_in_bin_5 + VLA_dI[5]
sim_data6 = I_in_bin_6 + VLA_dI[6]
sim_data7 = I_in_bin_7 + VLA_dI[7]
sim_data8 = I_in_bin_8 + VLA_dI[8]
sim_data9 = I_in_bin_9 + VLA_dI[9]
sim_data10 = I_in_bin_10 + VLA_dI[10]
sim_data11 = I_in_bin_11 + VLA_dI[11]
sim_data12 = I_in_bin_12 + VLA_dI[12]
sim_data13 = I_in_bin_13 + VLA_dI[13]
sim_data14 = I_in_bin_14 + VLA_dI[14]'''

sim_data = np.concatenate((sim_data3,sim_data4,sim_data5,sim_data6,sim_data7,sim_data8,sim_data9,sim_data10,sim_data11,sim_data12,sim_data13,sim_data14))
sim_data = sim_data[:,0]

plt.figure()
#plt.plot(M31_fr_pos[index_b2:index_b14],175*sim_data,'--')
#plt.errorbar(M31_fr_pos[index_b2::], 175*sim_data, yerr=175*randerr, fmt=".k", capsize=0)
plt.plot(M31_fr_pos,175*M31_ft_pos, 'r')
plt.xscale('log')
#plt.yscale('log')

plt.savefig('./sersic_profile_sim_errors.png', dpi=300, bbox_inches='tight')

##----------------------Guassian liklihood----------------------------------------------##

def log_likelihood(theta, x, y, yerr):
    m, b, log_f = theta
    model = m * x + b
    sigma2 = yerr**2 + model**2 * np.exp(2 * log_f)
    return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(sigma2))
