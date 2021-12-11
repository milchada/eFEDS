import numpy as np
import matplotlib.pylab as plt
from matplotlib import colors, cm
from astropy.io import fits
from astropy import units as u 
from astropy.cosmology import default_cosmology
lcdm = default_cosmology.get()

def make_arrays(prefix, energy_range='_0520',nbins=10, galsperfile=66): 
	#'_0510', '' = '_0520', '1020' = '0520' - '0510'
	
	files = glob.glob(prefix+'*%s*fits' % energy_range) 
	files.sort()
	cts = np.zeros((int(len(files)*galsperfile), nbins))
	exp = np.zeros((int(len(files)*galsperfile), nbins))
	area = np.zeros((int(len(files)*galsperfile), nbins))
	for i in range(len(files)):
		# try:
			f = fits.getdata(files[i])
			counts = f['COUNTS']
			length = int(len(counts)/nbins)
			shape = (length, nbins)
			cts[i*length:(i+1)*length] = counts.reshape(shape)
			exp[i*length:(i+1)*length] = f['MEAN_SRC_EXP'].reshape(shape)
			area[i*length:(i+1)*length] = (f['AREA'] / (20**2)).reshape(shape) #to convert pixels to arcsec
			print(i)
		# except (ValueError, KeyError):
		# 	print(i, ' oops!')
		# 	continue
	return cts, exp, area 

def format_plot(ax, yconv=3.1e39):
	ax.set_xlabel('R (kpc)', fontsize=14)
	ax.set_xlim(4,100)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_ylim(3e33, 4e37)
	ax.set_ylabel(r'erg/s/kpc$^2$', fontsize=14)
	plt.tight_layout()

from psf_contam import psf_curve, fit, xrb 
from psf_contam import area as area_psf
	
def profiles(cts, exp, area, sdss, yconv=3.1e39, bck=None, simsfile = '../Ms10.20_10.70.out', 
	ret_im = False, ret=True, bckcols=5, ax = None, xrb=False,color='tab:blue',cts_to_erg = 1.16e-12):
	#cts/s to erg/s/cm**2 default is for 0.5-2.0 keV input
	#if we're only using T5 and T7, we only use photons at 0.8-2.0keV, but we still want Lx in 0.5-2.0, 
	#so instead the conversion factor is 2.26e-12.
	
	z = sdss['z']
	dcorr = 4*np.pi*lcdm.luminosity_distance(z).to('cm')**2 #for luminosity distance
	dA = lcdm.angular_diameter_distance(z).to('kpc').value*u.arcsec.to('radian')
	
	area_kpc = np.zeros(area.shape)
	x = np.arange(0, cts.shape[1], 1)*10 + 5

	net_cts = np.zeros(cts.shape)
	for i in range(exp.shape[0]):
		area_kpc[i] = area[i]*(dA[i]**2)  #arcsec to kpc
		for j in range(exp.shape[1]):
			# bkg = bck['COUNTS'][i]/(bck['AREA'][i]/400.) * area[i,j] 
			bkg = np.nanmean(cts[i,-bckcols:]/area[i,-bckcols:]) * area[i,j]
			if not exp[i,j]:
				exp[i,j] = np.nan
				cts[i,j] = np.nan
			net_cts[i,j] = cts[i,j] - bkg

	sb = np.nansum(net_cts/area_kpc,axis=0) / np.nansum(exp, axis=0) 
	
	err = np.sqrt(np.nansum(cts, axis=0)) / np.nansum(net_cts, axis=0)
	errtot = np.sqrt(np.nansum(cts))/np.nansum((net_cts))

	#EVERYTHING SO FAR IS IN COUNTS
	sb *= cts_to_erg * np.mean(dcorr.value)
	if bck:
		sb_bkg = np.nansum(bck['COUNTS']/(bck['AREA']*dA**2 / 400.))
		sb_bkg /= (np.nansum(bck['MEAN_SRC_EXP']))
	else:
		sb_bkg = np.nanmean(cts[:,-bckcols:]/(area[:,-bckcols:]*exp[:,-bckcols:]))
	sb_bkg *= cts_to_erg * np.mean(dcorr.value)
	
	#XRB contribution
	if xrb:
		norm = sb_xrb(sdss).mean()/area_psf[0]
		yp = psf_curve(x, *fit)
		sbxrb = yp*norm
		sb -= sbxrb

	if ret_im:
		if not ax:
			fig, ax = plt.subplots()
		ax.errorbar(x, sb, xerr=5, yerr=err*sb, linewidth=0, elinewidth=1, color=color)
		
		if simsfile:
			ax.hlines(sb_bkg, 0,200, color='k', alpha=0.25)
			ax.hlines(0.05*sb_bkg, 0,200, color='k', alpha=0.25, linestyle='dotted')
			sims = np.genfromtxt(simsfile,skip_footer=27)
			ax.fill_between(sims[:,0],sims[:,2]*yconv,sims[:,3]*yconv,color='tab:blue',alpha=0.5)
			ax.fill_between(sims[:,0],sims[:,5]*yconv,sims[:,6]*yconv,color='tab:red',alpha=0.5)
			ax.fill_between(sims[:,0],sims[:,8]*yconv,sims[:,9]*yconv,color='tab:purple',alpha=0.5)
			ax.fill_between(sims[:,0],sims[:,-2]*yconv,sims[:,-1]*yconv,color='tab:orange',alpha=0.5)
		if xrb:
			x = np.arange(0,100,.1)
			y = psf_curve(x, *fit)
			print(x.shape, y.shape, norm.shape)
			ax.plot(x, y*norm, c='k', label=r'XRB $\times$ PSF')
		format_plot(ax)
		return ax 
	elif ret:
		return x, sb, errtot, sb_bkg	

