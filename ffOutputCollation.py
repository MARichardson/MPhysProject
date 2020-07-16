'''

This is an output collation code to be used in conjunction with firefly_marv.py
This written by Molly Richardson <UP824687__at__myport.ac.uk> as part of an
MPhys project with Daniel Thomas.

'''

import glob
import sys
import os
import os.path
import numpy as np
from astropy.io import fits

#allows arrays to be printed fully
np.set_printoptions(threshold=sys.maxsize)

#-------------------------------------------------------------------------------
# Enter directory path to firefly output folder
#-------------------------------------------------------------------------------

dirpath = '/Users/Molly_Richardson/Documents/UNIVERSITY/Physics Y4/MPhys Project/sciamaOut/8078-12703/'
plateIFU = os.path.basename(dirpath[:-1])

print(plateIFU)

#-------------------------------------------------------------------------------
# initialising arrays using IFU bundle to set size
#-------------------------------------------------------------------------------

temp = plateIFU.split('-')
plate = temp[0]
ifu = temp[1]
#print(ifu)

ifuDict = {'19': 34,
           '37': 44,
           '61': 54,
           '91': 64,
           '127': 74}

for key in ifuDict.keys():
    if key in ifu:
        size = ifuDict[key]
        print(size)

realsize = size -2
print(realsize)

array = (realsize,realsize)

lwmetal = np.zeros(array)
lwage = np.zeros(array)

#-------------------------------------------------------------------------------
# Extracting data from fits files
# Replace all instances of the plateIFU in the strings with your plateIFU
#-------------------------------------------------------------------------------

#Before extracting the data from the fits file, the coordinate of the file must
#be extrapolated from the file name, so the values for lwage and lwmetallicity
#can be appended into the arrays in the correct place. The technique used means
#the files can be imported in the incorrect order.

for filename in glob.glob(os.path.join(dirpath, '8078-12703_*.fits')):
    with fits.open(filename) as hdu:
        print(filename)
        y = filename.replace(dirpath,'') #removing path from filename e.g. 8078-12703_9-20.fits
        #print(y)
        s = os.path.splitext(y)[0] #removing .fits e.g. 9000-1901_9-20
        #print(s)
        b_string = s.replace('8078-12703_','') #removing plateIFU from name e.g. 9-20
        #print(b_string)
        a_string = b_string.replace('-',' ') #removing hyphen leaving coord with space as string e.g. 9 20
        #print(a_string)
        nums = [int(n) for n in a_string.split()] #turning string into array, elements seperated by spaces e.g. [9, 20]
        #print(nums)
        age = hdu[1].header['HIERARCH age_lightW'] #setting age variable as the age in the FIREFLY output fits file
        #print(age)
        lwage[nums[0],nums[1]] = age #appending the age to the correct point in the array defined above e.g. at [9, 20]
        metal = hdu[1].header['HIERARCH metallicity_lightW'] #repeating above process with metallicity
        #print(metal)
        lwmetal[nums[0],nums[1]] = metal
        hdu.close() #closing the fits file

#print(lwage)
#print(lwmetal)

#-------------------------------------------------------------------------------
# 2D density plots
#-------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np

mask2 = np.where(lwage == -9999) #removing error points which can be found in some galaxies
lwage[mask2] = 0
print(np.min(lwage))
plt.imshow(lwage, cmap='magma')
plt.title('Light-Weighted Age')
plt.xlabel('Spaxel')
plt.ylabel('Spaxel')
plt.gca().invert_yaxis()
plt.colorbar(label='Age log (age/Gyr)')
#Contours can be added by uncommenting the lines below
#x = np.arange(72)
#y = np.arange(72)
#plt.contour(x,y,lwage, cmap='GnBu')
#plt.colorbar(label='Age log (age/Gyr)')

import matplotlib.pyplot as plt
import numpy as np

#Same repeated below for metallicity map:

mask2 = np.where(lwmetal == -9999)
lwage[mask2] = 0
print(np.min(lwmetal))
plt.imshow(lwmetal, cmap='magma_r')
plt.title('Light-Weighted Metallicity')
plt.xlabel('Spaxel')
plt.ylabel('Spaxel')
plt.gca().invert_yaxis()
plt.colorbar(label='Metallcity [Z/H]')
#Contours can be added by uncommenting the lines below
#x = np.arange(72)
#y = np.arange(72)
#plt.contour(x,y,lwage, cmap='GnBu')
#plt.colorbar(label='Age log (age/Gyr)')

#-------------------------------------------------------------------------------
# lwage against radius plot
#-------------------------------------------------------------------------------

#Input dir path to maps file:
maps = fits.open('/Users/Molly_Richardson/Documents/UNIVERSITY/Physics Y4/MPhys Project/manga-8078-12703-MAPS-VOR10-MILESHC-MASTARHC.fits')

radius = maps[8].data
rad = radius[1,:,:]
#print(np.shape(rad))
#print(np.shape(lwage))

plt.scatter(rad,lwage, c='deepskyblue')
plt.plot([0,1.5],[0.6,0.59], c='purple',lw=3) #input LoBF grad here
grad1 = (0.59-0.6/1.5-0)
grad2 = format(grad1, '.6f')
gradient = 'gradient = {}'.format(grad2)
print(gradient)
plt.annotate(gradient, (1.5, 1))
plt.title('Light-Weighted Age')
plt.xlabel('Radius R/R$_{e}$')
plt.ylabel('log Age (Gyr)')
#print(rad[,12])

#Same repeated below for metallicity plot

plt.scatter(rad,lwmetal, c = 'deepskyblue')
plt.plot([0,1.5],[0.3,-0.6], c='purple',lw=3)
grad1 = (-0.6-0.3/2.2-0)
grad2 = format(grad1, '.6f')
gradient = 'gradient = {}'.format(grad2)
print(gradient)
plt.annotate(gradient, (1.5, 0.2))
plt.title('Light-Weighted Metallicity')
plt.xlabel('Radius R/R$_{e}$')
plt.ylabel('Metallcity [Z/H]')
