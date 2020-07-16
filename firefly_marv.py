
"""
.. moduleauthor:: Daniel Thomas <daniel.thomas__at__port.ac.uk>
.. contributions:: Johan Comparat <johan.comparat__at__gmail.com>
.. contributions:: Violeta Gonzalez-Perez <violegp__at__gmail.com>

Firefly is initiated with this script.
All input data and parmeters are now specified in this one file.

This wrapper script has been reformatted and now uses Marvin implementation,
this was done by Molly Richardson <UP824687__at__myport.ac.uk> as part of an
MPhys project with Daniel Thomas.

"""

import numpy as np
import sys, os
from os.path import join
import time
import firefly_setup as setup
import firefly_models as spm
import astropy.cosmology as co
import subprocess
import marvin
from marvin.tools import Maps
from marvin.tools import ModelCube

#-------------------------------------------------------------------------------
# To run this scipt:
# Input your plateIFU below (LINE 43)
# Input your file path to the instrumental resolution file LINE (112)
#
# Results will be output in a folder named after your plateIFU in your FIREFLY
# directory, each fits file will follow the naming convention
# plateIFU_xcoord-ycoord.fits
#-------------------------------------------------------------------------------

#allows arrays to be printed fully
np.set_printoptions(threshold=sys.maxsize)

#initiating time
t0=time.time()

input_file = '8078-12701' #input your plateIFU

#Load in data
maps = Maps(input_file, bintype='VOR10') #specifying Voronoi-cell binned data
mc = ModelCube(input_file, bintype='VOR10')
print('files found')

#RA and DEC
ra=float(maps.header['OBJRA']); dec=float(maps.header['OBJDEC'])

#Importing flux, emission, wavelength and error
findflux = mc.binned_flux
print('flux found')
findwavelength = mc.binned_flux.wavelength
wavelength = np.array(findwavelength)
print('wavelength found')
finderror = mc.binned_flux.ivar
print('error found')
findemission = mc.emline_fit
print('emission found')

#-------------------------------------------------------------------------------
# Naming convention for output files
#-------------------------------------------------------------------------------

#Provides spaxel coordinate to be used in output file name.

temp = input_file.split('-')
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

k=0

print(coord[k])


#-------------------------------------------------------------------------------
# Initiating variables and setting parameters
#-------------------------------------------------------------------------------

#Instrumental resolution
r_instrument = np.loadtxt(fname = "/Users/Molly_Richardson/Documents/UNIVERSITY/Physics Y4/MPhys Project/Firefly/firefly_release-0.1.1/MaNGA_spectral_resolution.txt")

#Masking emission lines
N_angstrom_masked=20 #N_angstrom_masked set to 20 in _init_ function

#Key which models and minimum age and metallicity of models to be used
models_key='m11'
cosmo = co.Planck15
ageMin = 0. ; #age max in loop
ZMin = 0.001 ; ZMax = 10.

#Model library3
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
#print(coord)

#setting i, j
size = len(maps['stellar_sigmacorr'])
#print(size)

for i in range(size):
#for i in np.arange(14, 73, 1):
    for j in range(size):
    #for j in np.arange(43, 74, 1):
        flux1 = np.array(findflux[:, i,j]) #flux
        #print(flux1)
        print(i)
        print(j)

        if np.max(flux1) > 0: #if flux value is not zero, produce output file

            emission = np.array(findemission[:, i,j])
            flux = flux1 - emission #flux-emission (possibly change names later on)
            error = 1/np.sqrt(finderror[:,i,j]) #1/sqrt because inverse variance
            #print(error)

            vel = maps['stellar_vel'].data[i,j] #finding velocity for redshift calulation
            redshift = maps.dapall['z'] + (vel/(2.99*10**8)) #redshift calculation
            #print(maps.dapall['z'])
            #print(redshift)

            restframe_wavelength = wavelength/(1+redshift)
            lines_mask = ((restframe_wavelength > 3728 - N_angstrom_masked) & (restframe_wavelength < 3728 + N_angstrom_masked)) | ((restframe_wavelength > 5007 - N_angstrom_masked) & (restframe_wavelength < 5007 + N_angstrom_masked)) | ((restframe_wavelength > 4861 - N_angstrom_masked) & (restframe_wavelength < 4861 + N_angstrom_masked)) | ((restframe_wavelength > 6564 - N_angstrom_masked) & (restframe_wavelength < 6564 + N_angstrom_masked))

            #Velocity dispersion in km/s
            vdisp = maps['stellar_sigma'].data[i,j]

            #Brought in from outside loop for redshift value
            ageMax = cosmo.age(redshift).value

            #Setting output file
            outputFolder = join(os.environ['FF_DIR'], input_file)
            #Setting file name using naming convention from outside loop
            output_file = join( outputFolder , input_file + '_' + coord[k] + ".fits")
            print(i)
            print(j)
            print(coord[k])
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
            spec=setup.firefly_setup(input_file,N_angstrom_masked=N_angstrom_masked)
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

        k=k+1 #Increasing k, naming convention variable, for next spaxel coordinate

print()
print ("Done... total time:", int(time.time()-t0) ,"seconds.")
print()
