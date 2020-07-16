'''
This code is designed to work with the results provided by FIREFLY.
This code produces a scatter plot of the data.
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


plateIFU = '8158-12702'

filename1 = join('/Users/Molly_Richardson/Documents/UNIVERSITY/Physics Y4/MPhys Project/finalResults/radialPlots3/' + plateIFU + 'lwagePlot3.png')
filename2 = join('/Users/Molly_Richardson/Documents/UNIVERSITY/Physics Y4/MPhys Project/finalResults/radialPlots3/' + plateIFU + 'lwmetalPlot3.png')

#-------------------------------------------------------------------------------
# Enter directory path to firefly output folder
#-------------------------------------------------------------------------------

dirpath = join('/Users/Molly_Richardson/Documents/UNIVERSITY/Physics Y4/MPhys Project/Firefly/firefly_release-0.1.1/Results/' + plateIFU + '/')
plateIFU = os.path.basename(dirpath[:-1])

print(plateIFU)

maps = fits.open(join('/Users/Molly_Richardson/Documents/UNIVERSITY/Physics Y4/MPhys Project/mapFiles/manga-' + plateIFU + '-MAPS-VOR10-MILESHC-MASTARHC.fits'))

radius = maps[8].data
rad = radius[1,:,:]

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

size = np.shape(rad)[0]

array = (size,size)

lwmetal = np.zeros(array)
lwage = np.zeros(array)


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
        removefits = os.path.splitext(removepath)[0] #removing .fits e.g. 9000-1901_9-20
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
# lwage against radius plot
#-------------------------------------------------------------------------------

plt.scatter(rad,lwage, c='tomato')
plt.plot([0,1.75],[0.72,0.57], c='indigo',lw=3) #input LoBF grad here
grad1 = ((0.57-0.72)/(1.75-0))
grad2 = format(grad1, '.6f')
gradient = 'gradient = {}'.format(grad2)
print(gradient)
plt.annotate(gradient, (0, -0.4))
plt.title(join(plateIFU + ' Light-Weighted Age'))
plt.xlabel('Radius R/R$_{e}$')
plt.ylabel('log Age (Gyr)')
plt.savefig(filename1)
plt.close()
#print(rad[,12])

#Same repeated below for metallicity plot

plt.scatter(rad,lwmetal, c = 'indigo')
plt.plot([0,1.75],[0.2,-0.73], c='tomato',lw=3)
grad1 = ((-0.73-0.2)/(1.75-0)) #y2-y1/x2-x1
grad2 = format(grad1, '.6f')
gradient = 'gradient = {}'.format(grad2)
print(gradient)
plt.annotate(gradient, (0, -1.7))
plt.title(join(plateIFU + ' Light-Weighted Metallicity'))
plt.xlabel('Radius R/R$_{e}$')
plt.ylabel('Metallcity [Z/H]')
plt.savefig(filename2)
