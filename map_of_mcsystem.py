#!/usr/bin/env/python

## Program to plot the Cepheids distances as a fn of position

import sys
sys.path.append('/Users/vs522/Dropbox/Python/')


import numpy as np
import matplotlib.pyplot as mp
import glob
import re
import os
from scipy.optimize import curve_fit
from matplotlib import cm
from scipy.interpolate import griddata
import numpy.ma as ma
import matplotlib.gridspec as gridspec


import coordinate_conversion as cc

from mpl_toolkits.mplot3d import Axes3D


os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin'
matplotlib.rc('text',usetex=True)
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Garamond']

files = glob.glob('*.glo_avs')


band = []
period = []
mag = []
err = []
cepid = []
ra = []
dec = []
cepname = []

## Coordintes of SMC centre
## A = 00:52:44.8
## D = -72:49:43

A, D = cc.todeg("00:52:44.8", "-72:49:43")

print(A, D)

Arad = A * np.pi / 180.
Drad = D * np.pi / 180.

for name in files:
	cepid.append(re.sub('.glo_avs', '', name))
	#print cepid
	for line in open(name,'r'):
		data1 = line.split()
		if data1[3] != '=' and data1[3] != 'std':
		## Need to remove the extra characters with regex -- re module
			band = (re.sub('[\[\]<>]', '', data1[3]))
			if band =='3.6':
				period.append(float(data1[2]))
				mag.append(float(data1[5]))
				if data1[7] != 'single':
					err.append(float(data1[9]))
				else:
					err.append(0.1)
				


logp = log10(period)
deltara = []
deltadec = []
radeg1 = []
decdeg1 = []

for coords in open("sorted_cepheid_coords","r"):
	data = coords.split()
	cepname.append(data[0])
	ra = data[1]+":"+data[2]+":"+data[3]
	dec = data[4]+":"+data[5]+":"+data[6]
	
	radeg, decdeg = cc.todeg(ra,dec)
	rarad = radeg * np.pi / 180.
	decrad = decdeg * np.pi / 180.
	radeg1.append(radeg)
	decdeg1.append(decdeg)
	
	deltara.append((radeg - A) * (cos(radians(decdeg))))
	deltadec.append(decdeg - D)
## Now do the LMC calculations

for lmccep in open("lmc_cepheids_data", "r"):
	data = lmccep.split()
	cepname.append(data[0])
	cepid.append(data[0])
	period.append(float(data[1]))
	mag.append(float(data[2]))
	err.append(float(data[3]))
	ra = data[4]
	dec = data[5]
	
	radeg, decdeg = cc.todeg(ra,dec)
	rarad = radeg * np.pi / 180.
	decrad = decdeg * np.pi / 180.
	radeg1.append(radeg)
	decdeg1.append(decdeg)
	
	deltara.append((radeg - A) * (cos(radians(decdeg))))
	deltadec.append(decdeg - D)


deltara = np.array(deltara)
deltadec = np.array(deltadec)	
cepname = np.array(cepname)
radeg1 = np.array(radeg1)
decdeg1 = np.array(decdeg1)

period = np.array(period)
mag = np.array(mag)
err = np.array(err)
cepid = np.array(cepid)

output = open("all_ceps_data", 'w')
output2 = open("all_ceps_true_pos","w")

distance = np.zeros(cepname.size)
disterr = np.zeros(cepname.size)
distmod_high = 0
distmod_low = 100
for cepheids in range (0, cepname.size):
	plval = -3.306 * (np.log10(period[cepheids]) - 1.0) -5.80
	# difference = measured - expected
	# positive = too faint = further away
	# negative = too bright = closer
	
	distmod = (mag[cepheids] - plval) 
	## Find max and minimum values
	if distmod > distmod_high:
		distmod_high = distmod
	if distmod < distmod_low:
		distmod_low = distmod
		
	point = np.where(cepname == cepid[cepheids])
	distance[point] = 10**((distmod + 5.)/5.) / 1000.
	disterr[point] = 0.461*distance[point]*sqrt(err[cepheids]**2 + 0.108**2)
	
	print("{0} {1} {2} {3} {4} {5} {6}".format(cepname[cepheids], deltara[cepheids], deltadec[cepheids], period[cepheids], mag[cepheids], err[cepheids], distance[cepheids]), file=output)
	print ("{0} {1} {2} {3} {4} {5} {6}".format(cepname[cepheids], radeg1[cepheids], decdeg1[cepheids], period[cepheids], mag[cepheids], err[cepheids], distance[cepheids]), file=output2)
	

output.close()
output2.close()
	





mp.clf()




fig = mp.figure()


ax3d = fig.add_subplot(111, projection="3d")

ax3d.scatter(deltara[(period < 60.) & (period > 6.)], distance[(period < 60.) & (period > 6.)], deltadec[(period < 60.) & (period > 6.)], c = distance[(period < 60.) & (period > 6.)],cmap=cm.Spectral_r, s=40)
ax3d.set_xlabel("$\Delta \\alpha$ (deg)")
ax3d.set_zlabel("$\Delta \\delta$ (deg)")
ax3d.set_ylabel("Distance (kpc)")


mp.show()

	
## fitting the tilt
	
