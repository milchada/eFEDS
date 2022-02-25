#so now the issue is that the point sources may be contaminating the signal
from astropy.io import fits
import numpy as np
import matplotlib.pylab as plt
from scipy.optimize import curve_fit
from scipy.integrate import quad

path='' #/data/gnarming/uchadaya/eFEDS/ if on hea comp
cdf = np.mean(np.mean(fits.getdata(path+'psf1.fits'), axis=1), axis=1)
eef = np.arange(0.4, 1.0, 0.05)
sb = eef/(np.pi*cdf**2)
area = np.pi * (cdf[1:]**2 - cdf[:-1]**2)
dE = eef[1:] - eef[:-1]
sb[1:] = dE / area
sb /= sb[0]
x = (cdf[1:] + cdf[:-1])/2
x = np.insert(x,0,cdf[0]/2)
x *= (4 * 1.45)

def psf_curve(r, A, rs, beta):
	return A * (1 + (r/rs)**2)**(-beta) #Moffat function, beta=2.0 fitted for XMM

def psf_area(r, A, rs, beta): 
	#integrating this from r1 to r2 should give me the fraction of photons enclosed between those radii
		#IF THE PSF IS PROPERLY NORMALISED
		#in other words, integrating this gives the factor that the raw PSF should be normalised by. 
	return 2*np.pi*r * psf_curve(r, A, rs, beta) 

fit, _ =  curve_fit(psf_curve, x, sb)
area_psf, _ = quad(psf_area, args=(fit[0],fit[1], fit[2]), a=0, b=np.infty) 
xp = np.arange(0,100,.1)
yp = psf_curve(xp, *fit)/area_psf

def Lxrb(sdss, alpha = 9.05e28, beta = 1.62e39, low=False, high=False, ism=False):  #alpha in erg/s/Msun, beta in erg/s/Msun/yr)
	# "Combine the outer few bins"
	# "See galaxies that host AGN based on eROSITA catalog"
	# #i.e. you might be stacking some faint AGN
	dalpha = 0.37e28
	dbeta = 0.22e39
	if low:
		alpha -= dalpha
		beta -= dbeta
	if high:
		alpha += dalpha 
		beta += dbeta
	if ism:
		beta += 8.3e38
	mstar = 10**sdss['logMass']
	sfr = 10**(sdss['ssfr'] - 9) * mstar #the -9 converts from Msun/Gyr to Msun/yr
	return alpha*mstar + beta*sfr

