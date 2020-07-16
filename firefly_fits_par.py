
"""
.. moduleauthor:: Daniel Thomas <daniel.thomas__at__port.ac.uk>
.. contributions:: Johan Comparat <johan.comparat__at__gmail.com>
.. contributions:: Violeta Gonzalez-Perez <violegp__at__gmail.com>

Firefly is initiated with this script.
All input data and parmeters are now specified in this one file.

This wrapper script has been reformatted and runs in full parallelisation, it is
designed to be used with fireflypar.sh and was originally run on the Sciama
supercomputer. This code was developed by Molly Richardson 
<UP824687__at__myport.ac.uk> as part of an MPhys project with Daniel Thomas.

"""

import numpy as np
import sys, os
from os.path import join
import time
import firefly_setup as setup
import firefly_models as spm
import astropy.cosmology as co
import subprocess
from astropy.io import fits
import gzip

#initiating time
t0=time.time()

jobid = os.environ['SLURM_ARRAY_TASK_ID']

#-------------------------------------------------------------------------------
# Input your directory path to your input file
#-------------------------------------------------------------------------------

logcubegz = gzip.open('/mnt/lustre/MaNGA/MPL-9/DAP/VOR10-MILESHC-MASTARHC/8078/12701/manga-7443-12703-LOGCUBE-VOR10-MILESHC-MASTARHC.fits.gz', 'rb')
logcube = fits.open(logcubegz)

mapsgz = gzip.open('/mnt/lustre/MaNGA/MPL-9/DAP/VOR10-MILESHC-MASTARHC/8078/12701/manga-7443-12703-MAPS-VOR10-MILESHC-MASTARHC.fits.gz', 'rb')
maps = fits.open(mapsgz)

dapall = fits.open('/mnt/lustre/MaNGA/MPL-9/dapall-v2_7_1-2.4.1.fits')

#RA and DEC
ra = logcube[0].header['OBJRA']; dec=logcube[0].header['OBJDEC']

#-------------------------------------------------------------------------------
# Naming convention for output files
#-------------------------------------------------------------------------------

#Provides spaxel coordinate to be used in output file name.

plateIFU = maps[0].header['PLATEIFU']
#print(plateIFU)


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

m = []
for i in range(size):
    m.append(str(i))
#print("m = %s" % m)
n = []
for i in range(size):
    n.append(str(i))
#print("n = %s" % n)

#print(m[0])
coord = []
for p in range(size):
    for q in range(size):
        coord.append(m[p] + '-' + n[q])

#print(coord)

#k=0
#print(k)
#print(coord[k])

#-------------------------------------------------------------------------------
# Initiating variables and setting parameters
#-------------------------------------------------------------------------------

#Instrumental resolution
r_instrument = np.loadtxt(fname = "/users/mr824687/firefly_release-0.1.1/MaNGA_spectral_resolution.txt")

#Masking emission lines
N_angstrom_masked=20 #N_angstrom_masked set to 20 in _init_ function

#Key which models and minimum age and metallicity of models to be used
models_key='m11'
cosmo = co.Planck15
ageMin = 0. ; #age max in loop
ZMin = 0.001 ; ZMax = 10.

#Model library
model_libs=['MILES']

#Model imf
imfs=['kr']

#Specify whether data in air or vaccum
data_wave_medium='vacuum'

#Firefly assumes flux units of erg/s/A/cm^2. Need to choose factor in case flux is scaled
flux_units=10**(-17) #flux_units=10**(-17) for SDSS

#specify whether models should be downgraded to the instrumental resolution and galaxy velocity dispersion
downgrade_models=True

#specify whether write results
write_results=True

#-------------------------------------------------------------------------------
# Nested Loops
#-------------------------------------------------------------------------------

#setting i, j
#size = len(maps[18].data)
#print(size)
size = len(logcube[1].data[0,:,:])
#print(size)

#print(size)

for j in range(size):

    flux1 = logcube[1].data[:,jobid,j] #flux
    #print(i)
    #print(j)
    #print(coord[k])

    if np.max(flux1) > 0: #if flux value is not zero, produce output file

        #Importing flux, emission, wavelength and error
        flux = logcube[1].data[:,jobid,j] - logcube[8].data[:,jobid,j]
        wavelength = logcube[4].data
        error = 1/np.sqrt(logcube[2].data[:,jobid,j]) #1/sqrt because inverse variance

        vel = maps[15].data[jobid,j] #finding velocity for redshift calulation
        index = np.where(dapall[1].data['PLATEIFU'] == '7443-12703') #finding redshift
        redshift = dapall[1].data['NSA_Z'][index][0] + (vel/(2.99*10**8)) #redshift calculation
        #print(maps.dapall['z'])
        #print(redshift)

        restframe_wavelength = wavelength/(1+redshift)
        lines_mask = ((restframe_wavelength > 3728 - N_angstrom_masked) & (restframe_wavelength < 3728 + N_angstrom_masked)) | ((restframe_wavelength > 5007 - N_angstrom_masked) & (restframe_wavelength < 5007 + N_angstrom_masked)) | ((restframe_wavelength > 4861 - N_angstrom_masked) & (restframe_wavelength < 4861 + N_angstrom_masked)) | ((restframe_wavelength > 6564 - N_angstrom_masked) & (restframe_wavelength < 6564 + N_angstrom_masked))

        #Velocity dispersion in km/s
        vdisp = maps[18].data[jobid,j]

        #Brought in from outside loop for redshift value
        ageMax = cosmo.age(redshift).value

        #Setting output file
        outputFolder = join( os.environ['FF_DIR'], plateIFU)
        #Setting file name using naming convention from outside loop
        output_file = join( outputFolder , plateIFU + '_' + str(jobid) + '-' + str(j) + '.fits')
        #print(i)
        #print(j)
        #print(coord[k])
        #Duplicate file check
        if os.path.isfile(output_file):
            print()
            print('Warning: This object has already been processed, the file will be over-witten.')
            answer = input('** Do you want to continue? (Y/N)')
            if (answer=='N' or answer=='n'):
                sys.exit()
            os.remove(output_file)
        #Creating output fits file
        if os.path.isdir(outputFolder)==False:
            os.mkdir(outputFolder)
        print()
        print( 'Output file: ', output_file                 )
        print()
        prihdr = spm.pyfits.Header()
        prihdr['FILE']      = os.path.basename(output_file)
        prihdr['MODELS']    = models_key
        prihdr['FITTER']    = "FIREFLY"
        prihdr['AGEMIN']    = str(ageMin)
        prihdr['AGEMAX']    = str(ageMax)
        prihdr['ZMIN']      = str(ZMin)
        prihdr['ZMAX']      = str(ZMax)
        prihdr['redshift']  = redshift
        prihdr['HIERARCH age_universe']    = np.round(cosmo.age(redshift).value,3)
        prihdu = spm.pyfits.PrimaryHDU(header=prihdr)
        tables = [prihdu]

        #define input object to pass data on to firefly modules and initiate run
        spec=setup.firefly_setup(plateIFU,N_angstrom_masked=N_angstrom_masked)
        spec.openSingleSpectrum(wavelength, flux, error, redshift, ra, dec, vdisp, lines_mask, r_instrument)
        #spec.openMANGASpectrum(data_release, path_to_logcube, path_to_drpall, bin_number, plate_number, ifu_number)

        #Error handling - did not converge = could not produce fit
        did_not_converge = 0.
        #print(did_not_converge)
        try :
            #prepare model templates
            model = spm.StellarPopulationModel(spec, output_file, cosmo, models = models_key, model_libs = model_libs, imfs = imfs, age_limits = [ageMin,ageMax], downgrade_models = downgrade_models, data_wave_medium = data_wave_medium, Z_limits = [ZMin,ZMax], use_downgraded_models = False, write_results = write_results, flux_units=flux_units)
            #initiate fit
            model.fit_models_to_data()
            tables.append( model.tbhdu )
        except (ValueError):
            tables.append( model.create_dummy_hdu() )
            did_not_converge +=1
            print('did not converge')
        if did_not_converge < 1 :
            complete_hdus = spm.pyfits.HDUList(tables)
            if os.path.isfile(output_file):
                os.remove(output_file)
            complete_hdus.writeto(output_file)

    #k=k+1 #Increasing k, naming convention variable, for next spaxel coordinate

print()
print ("Done... total time:", int(time.time()-t0) ,"seconds.")
print()
