import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import astropy.units as u

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
zvals = np.arange(0, 6, 0.1)
dist = cosmo.angular_diameter_distance(zvals)

plt.rc('xtick.major', size=4)
plt.rc('ytick.major', size=4)
plt.rc('xtick.minor', size=2)
plt.rc('ytick.minor', size=2)
plt.rc('axes', grid=False)
plt.rc('xtick.major', width=1)
plt.rc('xtick.minor', width=1)
plt.rc('ytick.major', width=1)
plt.rc('ytick.minor', width=1)

from astropy.cosmology import z_at_value
ages = np.array([13, 10, 8, 6, 5, 4, 3, 2, 1.5, 1.2, 1])*u.Gyr
ageticks = [z_at_value(cosmo.age, age) for age in ages]

from astropy.cosmology import Planck13
dist2 = Planck13.angular_diameter_distance(zvals)

fig = plt.figure(figsize=(6,4))
ax = fig.add_subplot(111)
ax.plot(zvals, dist2, label='Planck 2013')
ax.plot(zvals, dist, label='$h=0.7,\ \Omega_M=0.3,\ \Omega_\Lambda=0.7$')
ax.legend(frameon=0, loc='lower right')
ax2 = ax.twiny()
ax2.set_xticks(ageticks)
ax2.set_xticklabels(['{:g}'.format(age) for age in ages.value])
zmin, zmax = 0.0, 5.9
ax.set_xlim(zmin, zmax)
ax2.set_xlim(zmin, zmax)
ax2.set_xlabel('Time since Big Bang (Gyr)')
ax.set_xlabel('Redshift')
ax.set_ylabel('Angular diameter distance (Mpc)')
ax.minorticks_on()
ax.set_ylim(0, 1890)
fig.savefig('ang_dist.png', dpi=200, bbox_inches='tight')
