'''
This code is designed to work with the results provided by FIREFLY.
This code produces a linear histogram plot of the data.
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


plateIFU = '8154-12704'

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
# lwage against radius plot from grad
#-------------------------------------------------------------------------------

rad2=rad.flatten()
lwage2=lwage.flatten()
lwmetal2=lwmetal.flatten()

mask3 = np.where(lwage2 > 1500) #removing error points which can be found in some galaxies
lwage2[mask3] = 0

plt.figure()
x=rad2
y=lwage2

grad=-0.0212333
err=0.049
x1=0.67
y1=0
y2=1.5
x2=((y2-y1)*grad)+x1


plt.xlabel('Radius R/R$_{e}$', fontsize=13)
plt.ylabel('log Age (Gyr)', fontsize=13)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.title(join(plateIFU+' Light-Weighted Age'), fontsize=14)
plt.scatter(x,y,s=0.75,color='k')
h =plt.hist2d(x, y, bins=25, cmin=5, cmax=300, cmap='plasma')
cb=plt.colorbar(h[3])
plt.plot([y1,y2],[x1,x2], c='k',lw=3) #input LoBF grad here
gradient = join(r'$\nabla$ Age = {}'.format(grad) + r'$\pm${}'.format(err))
agan=plt.annotate(gradient, (0.05, -0.34), fontsize=13)
plt.savefig(filename1)

#-------------------------------------------------------------------------------
# lwmetal against radius plot from grad
#-------------------------------------------------------------------------------
#
# rad2=rad.flatten()
# lwage2=lwage.flatten()
# lwmetal2=lwmetal.flatten()
#
# mask3 = np.where(lwage2 > 1500) #removing error points which can be found in some galaxies
# lwage2[mask3] = 0
#
# grad=-0.201001
# err=0.085
# x1=-0.03
# y1=0
# y2=1.5
# x2=((y2-y1)*grad)+x1
#
# plt.figure()
# x=rad2
# y=lwmetal2
# plt.xlabel('Radius R/R$_{e}$', fontsize=13)
# plt.ylabel('Metallcity [Z/H]', fontsize=13)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.title(join(plateIFU+' Light-Weighted Metallicity'), fontsize=14)
# plt.scatter(x,y,s=0.75,color='k')
# h =plt.hist2d(x, y, bins=25, cmin=10, cmax=300, cmap='plasma')
# cb=plt.colorbar(h[3])
# plt.plot([y1,y2],[x1,x2], c='k',lw=3) #input LoBF grad here
# gradient = join(r'$\nabla$ [Z/H] = {}'.format(grad) + r'$\pm${}'.format(err))
# agan=plt.annotate(gradient, (0.05, -1), fontsize=13)
# plt.savefig(filename2)

#-------------------------------------------------------------------------------
# lwage against radius plot
#-------------------------------------------------------------------------------
#
# rad2=rad.flatten()
# lwage2=lwage.flatten()
# lwmetal2=lwmetal.flatten()
# #np.savetxt("rad.csv", rad, delimiter=",")
# #np.savetxt("lwage.csv", lwage, delimiter=",")
#
# #write(rad2)
#
# mask3 = np.where(lwage2 > 1500) #removing error points which can be found in some galaxies
# lwage2[mask3] = 0
#
# plt.figure()
# x=rad2
# y=lwage2
# #plotrange=[-0.5, 5.5]
# plt.xlabel('Radius R/R$_{e}$', fontsize=13)
# plt.ylabel('log Age (Gyr)', fontsize=13)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.title(join(plateIFU+' Light-Weighted Age'), fontsize=14)
# plt.scatter(x,y,s=0.75,color='k')
# #plt.scatter(x,y,color='k')
# h =plt.hist2d(x, y, bins=25, cmin=5, cmax=300, cmap='plasma')
# #h =plt.hist2d(x, y, range=[plotrange,plotrange], bins=50)
# #h = plt.hist2d(x, y)
# cb=plt.colorbar(h[3])
# plt.plot([0,1.5],[0.4,0.2], c='k',lw=3) #input LoBF grad here
# grad1 = ((0.2-0.4)/(1.5-0))
# grad2 = format(grad1, '.6f')
# gradient = r'$\nabla$ Age = {}'.format(grad2)
# agan=plt.annotate(gradient, (0.05, -0.6), fontsize=13)
#
#
# #figure.set_size_inches(6, 4)
# #xl=np.arange(0,0.4,0.01,dtype=float)
# #plt.scatter(xl,xl,color='k',marker=".",s=2)
# #plt.savefig(filename1)
#
# #cb.remove()
# #agan.remove()
#
# plt.figure()
# x=rad2
# y=lwmetal2
# #plotrange=[-0.5, 5.5]
# plt.xlabel('Radius R/R$_{e}$', fontsize=13)
# plt.ylabel('Metallcity [Z/H]', fontsize=13)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.title(join(plateIFU+' Light-Weighted Metallicity'), fontsize=14)
# plt.scatter(x,y,s=0.75,color='k')
# #plt.scatter(x,y,color='k')
# h =plt.hist2d(x, y, bins=25, cmin=5, cmax=300, cmap='plasma')
# #h =plt.hist2d(x, y, range=[plotrange,plotrange], bins=50)
# #h = plt.hist2d(x, y)
# plt.colorbar(h[3])
# plt.plot([0,1.5],[0,-0.4], c='k',lw=3) #input LoBF grad here
# grad1 = ((-0.4-0)/(1.5-0))
# grad2 = format(grad1, '.6f')
# gradient = r'$\nabla$ [Z/H] = {}'.format(grad2)
# plt.annotate(gradient, (0.05,-0.5), fontsize=13)
# figure2 = plt.gcf()
# #figure.set_size_inches(6, 4)
# #plt.savefig(filename2)
#
# plt.show()
#
#
