'''
This code is designed to work with the results provided by FIREFLY.
This code produces a 2D map plot of the data.
Developed by Molly Richardson <UP824687__at__myport.ac.uk>
'''

import glob
import sys
import os
import os.path
from os.path import join
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np


plateIFU = input("Input desired plateIFU: ")

filename1 = join('/Users/Molly_Richardson/Documents/UNIVERSITY/Physics Y4/MPhys Project/finalResults/mapPlots3/' + plateIFU + 'lwageMap.png')
filename2 = join('/Users/Molly_Richardson/Documents/UNIVERSITY/Physics Y4/MPhys Project/finalResults/mapPlots3/' + plateIFU + 'lwageMapSN10.png')
filename3 = join('/Users/Molly_Richardson/Documents/UNIVERSITY/Physics Y4/MPhys Project/finalResults/mapPlots3/' + plateIFU + 'lwageMapSN20.png')
filename4 = join('/Users/Molly_Richardson/Documents/UNIVERSITY/Physics Y4/MPhys Project/finalResults/mapPlots3/' + plateIFU + 'lwmetalMap.png')
filename5 = join('/Users/Molly_Richardson/Documents/UNIVERSITY/Physics Y4/MPhys Project/finalResults/mapPlots3/' + plateIFU + 'lwmetalMapSN10.png')
filename6 = join('/Users/Molly_Richardson/Documents/UNIVERSITY/Physics Y4/MPhys Project/finalResults/mapPlots3/' + plateIFU + 'lwmetalMapSN20.png')

#-------------------------------------------------------------------------------
# Enter directory path to firefly output folder
#-------------------------------------------------------------------------------

dirpath = join('/Users/Molly_Richardson/Documents/UNIVERSITY/Physics Y4/MPhys Project/Firefly/firefly_release-0.1.1/Results/' + plateIFU + '/')
plateIFU = os.path.basename(dirpath[:-1])

#dirpath = join('/Users/Molly_Richardson/Documents/UNIVERSITY/Physics Y4/MPhys Project/Firefly/firefly_release-0.1.1/Results/' + plateIFU + '_mask/')
#plateIFU = os.path.basename(dirpath[:-6])


print(plateIFU)

maps = fits.open(join('/Users/Molly_Richardson/Documents/UNIVERSITY/Physics Y4/MPhys Project/mapFiles/manga-' + plateIFU + '-MAPS-VOR10-MILESHC-MASTARHC.fits'))
#maps.info()
#maps[1].header

#HALPHA

halpha = maps[30].data[24,:,:]

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
        #print(size)


size = np.shape(halpha)[0]

#array = (realsize,realsize)
array = (size,size)


lwmetal = np.zeros(array)
lwage = np.zeros(array)

#np.shape(lwage)

#-------------------------------------------------------------------------------
# Extracting data from fits files
#-------------------------------------------------------------------------------

#Before extracting the data from the fits file, the coordinate of the file must
#be extrapolated from the file name, so the values for lwage and lwmetallicity
#can be appended into the arrays in the correct place. The technique used means
#the files can be imported in the incorrect order.

for filename in glob.glob(os.path.join(dirpath, join(plateIFU + '_*.fits'))):
    with fits.open(filename) as hdu:
        #print(filename)
        removepath = filename.replace(dirpath,'') #removing path from filename e.g. 8078-12703_9-20.fits
        #print(removepath)

        removemask = removepath.replace('_mask','') #removing .fits e.g. 9000-1901_9-20
        #print(removemask)
        removefits = os.path.splitext(removemask)[0] #removing .fits e.g. 9000-1901_9-20

        #removefits = os.path.splitext(removepath)[0] #removing .fits e.g. 9000-1901_9-20
        #print(removefits)
        removeplateifu = removefits.replace(join(plateIFU +'_'),'') #removing plateIFU from name e.g. 9-20
        #print(removeplateifu)
        replacehyphen = removeplateifu.replace('-',' ') #removing hyphen leaving coord with space as string e.g. 9 20
        #print(replacehyphen)
        arraypos = [int(n) for n in replacehyphen.split()] #turning string into array, elements seperated by spaces e.g. [9, 20]
        #print(arraypos)
        age = hdu[1].header['HIERARCH age_lightW'] #setting age variable as the age in the FIREFLY output fits file
        #print(age)
        lwage[arraypos[0],arraypos[1]] = age #appending the age to the correct point in the array defined above e.g. at [9, 20]
        metal = hdu[1].header['HIERARCH metallicity_lightW'] #repeating above process with metallicity
        #print(metal)
        lwmetal[arraypos[0],arraypos[1]] = metal
        hdu.close() #closing the fits file

#-------------------------------------------------------------------------------
# LIGHT WEIGHTED AGE
#-------------------------------------------------------------------------------

masknan=np.where(halpha==0)
lwage[masknan]=float('Nan')
#plt.imshow(lwage, cmap='magma', origin='lower')


#HALPHA MASKING

halpha[32,59]=0
mask2 = np.where(halpha > 45) #removing error points which can be found in some galaxies
halpha[mask2] = 0
mask3 = np.where(halpha < -600) #removing error points which can be found in some galaxies
halpha[mask3] = 0

#LWAGE MAP ---------------------------------------------------------------------

mask4 = np.where(lwage == -9999) #removing error points which can be found in some galaxies
lwage[mask4] = 0
#print(np.min(lwage))
plt.imshow(lwage, cmap='magma')
plt.title(join(plateIFU + ' Light-Weighted Age'))
plt.xlabel('Spaxel')
plt.ylabel('Spaxel')
plt.gca().invert_yaxis()
cb=plt.colorbar(label='Age log (age/Gyr)')
x = np.arange(size)
y = np.arange(size)
hal=plt.contour(x,y,halpha,  cmap='Greys_r')
cb2=plt.colorbar(label='H-alpha 6564 ($1 \cdot 10^{-17}$ergs$^{-1}$spaxel$^{-1}$cm$^{-2}$)')
plt.savefig(filename1)
cb.remove()
cb2.remove()
plt.clf()


#SETTING S/N >10

SN1= maps[14].data
SN = np.where(SN1 < 10)

lwageSN = lwage
lwageSN[SN] = 0


#LWAGE MAP S/N >10 -------------------------------------------------------------
maskzero=np.where(lwage==0)
lwage[maskzero]=float('Nan')
plt.imshow(lwage, cmap='magma')
plt.title(join(plateIFU + ' Light-Weighted Age (S/N >10)'))
plt.xlabel('Spaxel')
plt.ylabel('Spaxel')
plt.gca().invert_yaxis()
cb=plt.colorbar(label='Age log (age/Gyr)')
x = np.arange(size)
y = np.arange(size)
hal=plt.contour(x,y,halpha, cmap='Greys_r')
cb2=plt.colorbar(label='H-alpha 6564 ($1 \cdot 10^{-17}$ergs$^{-1}$spaxel$^{-1}$cm$^{-2}$)')
plt.savefig(filename2)
cb.remove()
cb2.remove()
plt.clf()

#SETTING S/N >20

SN1= maps[14].data
SN = np.where(SN1 < 20)

lwageSN = lwage
lwageSN[SN] = 0


#LWAGE MAP S/N >20 -------------------------------------------------------------
maskzero=np.where(lwage==0)
lwage[maskzero]=float('Nan')
plt.imshow(lwage, cmap='magma')
plt.title(join(plateIFU + ' Light-Weighted Age (S/N >20)'))
plt.xlabel('Spaxel')
plt.ylabel('Spaxel')
plt.gca().invert_yaxis()
cb=plt.colorbar(label='Age log (age/Gyr)')
x = np.arange(size)
y = np.arange(size)
hal=plt.contour(x,y,halpha, cmap='Greys_r')
cb2=plt.colorbar(label='H-alpha 6564 ($1 \cdot 10^{-17}$ergs$^{-1}$spaxel$^{-1}$cm$^{-2}$)')
plt.savefig(filename3)
cb.remove()
cb2.remove()


#-------------------------------------------------------------------------------
# LIGHT WEIGHTED METAL
#-------------------------------------------------------------------------------

#FINDING D4000
#
# D4000 = maps[49].data[44,:,:]
# #plt.imshow(D4000)
# #print(D4000)
# D4000[np.isnan(D4000)]=0
# #print(D4000)
# #plt.imshow(D4000)
# masknan=np.where(D4000==0)
# lwmetal[masknan]=float('Nan')
# #plt.imshow(lwmetal, cmap='magma_r', origin='lower')
#
#
# #plt.imshow(D4000, origin='lower')
#
# #MASKING D4000
#
# #LWMETAL MAP -------------------------------------------------------------------
#
# #mask2 = np.where(lwmetal == -9999)
# #lwage[mask2] = 0
# #print(np.min(lwmetal))
# plt.imshow(lwmetal, cmap='magma_r')
# plt.title(join(plateIFU + ' Light-Weighted Metallicity'))
# plt.xlabel('Spaxel')
# plt.ylabel('Spaxel')
# plt.gca().invert_yaxis()
# cb=plt.colorbar(label='Metallcity [Z/H]')
# x = np.arange(size)
# y = np.arange(size)
# plt.contour(x,y,D4000, cmap='Greys_r')
# cb2=plt.colorbar(label='D4000')
# plt.savefig(filename4)
# cb.remove()
# cb2.remove()
# plt.clf()
#
# #SETTING S/N >10
# #
# # SN1= maps[14].data
# # SN = np.where(SN1 < 10)
# #
# # lwmetalSN = lwmetal
# # lwmetalSN[SN] = 0
# #
# # #LWMETAL MAP S/N>10 ------------------------------------------------------------
# #
# # mask2 = np.where(lwmetal == -9999)
# # lwage[mask2] = 0
# # print(np.min(lwmetal))
# # plt.imshow(lwmetal, cmap='magma_r')
# # plt.title(join(plateIFU + ' LW Metallicity (S/N >10)'))
# # plt.xlabel('Spaxel')
# # plt.ylabel('Spaxel')
# # plt.gca().invert_yaxis()
# # cb=plt.colorbar(label='Metallcity [Z/H]')
# # x = np.arange(74)
# # y = np.arange(74)
# # plt.contour(x,y,D4000, [0.25,0.75,1.25,1.5,1.75], cmap='Greys_r')
# # cb2=plt.colorbar(label='D4000')
# # plt.savefig(filename5)
# # cb.remove()
# # cb2.remove()
# #
# # #SETTING S/N >20
# #
# # SN1= maps[14].data
# # SN = np.where(SN1 < 20)
# #
# # lwmetalSN = lwmetal
# # lwmetalSN[SN] = 0
# #
# # #LWMETAL MAP S/N>20 ------------------------------------------------------------
# #
# # mask2 = np.where(lwmetal == -9999)
# # lwage[mask2] = 0
# # print(np.min(lwmetal))
# # plt.imshow(lwmetal, cmap='magma_r')
# # plt.title(join(plateIFU + ' LW Metallicity (S/N >20)'))
# # plt.xlabel('Spaxel')
# # plt.ylabel('Spaxel')
# # plt.gca().invert_yaxis()
# # cb=plt.colorbar(label='Metallcity [Z/H]')
# # x = np.arange(74)
# # y = np.arange(74)
# # plt.contour(x,y,D4000, [0.25,0.75,1.25,1.5,1.75], cmap='Greys_r')
# # cb2=plt.colorbar(label='D4000')
# # plt.savefig(filename6)
# # cb.remove()
# # cb2.remove()
