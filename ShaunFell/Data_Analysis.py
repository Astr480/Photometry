import astropy
from astropy.io import fits
import numpy as np
import glob
import ccdproc
import photutils
from photutils import DAOStarFinder, CircularAperture, aperture_photometry, EllipticalAnnulus, CircularAnnulus
from photutils import psf; 
from photutils.psf import BasicPSFPhotometry
from photutils.psf.sandbox import DiscretePRF
from photutils.detection import IRAFStarFinder
from photutils.psf import IntegratedGaussianPRF, DAOGroup
from photutils.background import MMMBackground, MADStdBackgroundRMS
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.stats import SigmaClip, mad_std
from photutils import MedianBackground
import matplotlib.pyplot as plt
from matplotlib.colors import Colormap
from astropy.visualization import ZScaleInterval, ImageNormalize
import argparse

########### Argument parsing for command line arguments ########
P = argparse.ArgumentParser()
P.add_argument('-i', '--images', action = 'store_true', default=False)
args = P.parse_args()
show_images = args.images
################################################################


origin = 'lower'
cmp = 'gray'

dir = "D:/Users/shaun/Documents/01 UofW/Current Classes/Astro 480/Data Reduction/"
bias_files = glob.glob('Wolf1346/Bias*.fits')
wolf_files = glob.glob('Wolf1346/wolf*.fits')
flat_files = glob.glob('Wolf1346/flat*.fits')

biass = list(map(fits.open, bias_files))
wolf = list(map(fits.open, wolf_files))
flats = list(map(fits.open, flat_files))


bias_data = [np.asarray(i[0].data) for i in biass]
wolf_data = [np.asarray(i[0].data) for i in wolf]
flat_data = [np.asarray(i[0].data) for i in flats]

master_bias = sum(bias_data)/len(bias_data)

wolf_data_minus_bias = [i - master_bias for i in wolf_data]

if show_images:
	norm = ImageNormalize(wolf_data_minus_bias[0], interval=ZScaleInterval())
	plt.imshow(wolf_data_minus_bias[0], norm = norm, origin=origin, cmap = cmp); plt.show()
	plt.clf()



print("Averages of each of the bias frames: ", list(map(np.mean, bias_data)))
print("\n\nAverage of the master bias: ", np.mean(master_bias))


"""Skipping overscan
Xmask = np.arange(1026, 1077) #[1026:1077,]
Ymask = np.arange(1024, 1026) #[:,1024:1026]
wolf_data_minus_bias_minus_overscan = [np.delete(np.delete(i, Xmask, 0), Ymask, 1) for i in wolf_data_minus_bias]
if show_images:
	norm = ImageNormalize(wolf_data_minus_bias_minus_overscan[0], interval=ZScaleInterval())
	plt.imshow(wolf_data_minus_bias_minus_overscan[0], norm = norm, origin = origin, cmap = cmp); plt.show()
	plt.clf()
"""

master_flat = sum(flat_data)/len(flat_data)

print("I have the filter 'wash m' and 'J-C Rc")
for i in flats:
	if i[0].header['Filter'] == 'Wash M':
		print('The File {0} is using the wash m filter'.format(i[0].header['FILENAME']))
	elif i[0].header['Filter'] == 'J-C Rc':
		print('The File {0} is using the J-C Rc filter'.format(i[0].header['FILENAME']))








#Photometry

##identity stars
quadRU = fits.open('quadRU.fits')
bias = 100
fwhm = 5
quadRU_data = quadRU[0].data
if show_images:
	norm = ImageNormalize(quadRU_data, interval=ZScaleInterval())
	plt.imshow(quadRU_data, origin='lower', norm=norm, cmap='BrBG', clim=(0,1000)); plt.show()
	plt.clf()

DAO_stars = DAOStarFinder(bias, fwhm)
stars = DAO_stars.find_stars(quadRU_data)
print("Identified stars", stars)

##Use circular aperature
star_coords = [(stars['xcentroid'][i], stars['ycentroid'][i]) for i in range(len(stars['id']))]
apertures = CircularAperture(star_coords, r=3.)
print("Apertures", apertures)
phot_table = aperture_photometry(quadRU_data, apertures, method='exact')
print("Phot_table", phot_table)

if show_images:
	apertures.plot(color='blue', lw=2); 
	norm = ImageNormalize(quadRU_data, interval=ZScaleInterval())
	plt.imshow(quadRU_data, origin='lower', norm=norm, cmap='BrBG', clim=(0,1000)); plt.show()
	plt.clf()


##Use EllipticalAnnulus
annuli = EllipticalAnnulus(star_coords, a_in = 5.5, a_out = 7, b_out = 6, theta = np.pi/4)
print("Elliptical Annuli", annuli)
phot_table_ellip = aperture_photometry(quadRU_data, annuli, method='exact')
print("phot_table_ellip", phot_table_ellip)

if show_images:
	annuli.plot(color='blue', lw=2)
	plt.imshow(quadRU_data, origin = 'lower', norm = norm, cmap = 'BrBG', clim=(0,1000)); plt.show()
	plt.clf()
	
	
##use CircularAnnulus
annuli_c = CircularAnnulus(star_coords, r_in = 5.5, r_out = 7)
print("Circular Annuli", annuli_c)
phot_table_circ = aperture_photometry(quadRU_data, annuli_c, method='exact')
print("phot_table_circ", phot_table_circ)
if show_images:
	annuli_c.plot(color='blue', lw=2)
	plt.imshow(quadRU_data, origin='lower', norm = norm, cmap = 'BrBG', clim=(0,1000)); plt.show()
	plt.clf()

##background subtraction

##ellip annuli background subtraction
apers = [apertures, annuli, annuli_c]
phot_table = aperture_photometry(quadRU_data, apers)
background_mean_ellip = phot_table['aperture_sum_1']/annuli.area()
background_sum_ellip = background_mean_ellip*apertures.area()
final_sum_ellip = phot_table['aperture_sum_0'] - background_sum_ellip
phot_table['residual_aperture_sum'] = final_sum_ellip
print(phot_table['residual_aperture_sum'])

##circ annuli background subtraction
apers = [apertures, annuli, annuli_c]
phot_table = aperture_photometry(quadRU_data, apers)
background_mean_circ = phot_table['aperture_sum_2']/annuli.area()
background_sum_circ = background_mean_circ*apertures.area()
final_sum_circ = phot_table['aperture_sum_0'] = background_sum_circ
phot_table['residual_aperture_sum'] = final_sum_circ
print(phot_table['residual_aperture_sum'])


##background comparison
if show_images:
	plt.title("Background comparison between circular annuli and elliptic annuli")
	plt.ylabel("Residual aperture sum")
	plt.plot(range(len(final_sum_ellip)), final_sum_ellip, color='blue')
	plt.plot(range(len(final_sum_circ)), final_sum_circ, color='black')
	plt.show()
	plt.clf()



#PSF Photometry
"""I could not figure this part out before the due date
sigma =2.0
bkg = MedianBackground(sigma)
daogroup = DAOGroup(2.0*sigma*gaussian_sigma_to_fwhm)
mmm_bkg = MMMBackground()
fitter = LevMarLSQFitter()
psf_model = IntegratedGaussianPRF(sigma=sigma)
PSF_phot = BasicPSFPhotometry(group_maker=daogroup, bkg_estimator = mmm_bkg, psf_model = psf_model, fitter = LevMarLSQFitter(),fitshape = (11,11))
result_tab = PSF_phot(image=quadRU_data, init_guesses=stars)
residual_images = PSF_phot.get_residual_image()
plt.imshow(residual_images); plt.show()
"""





