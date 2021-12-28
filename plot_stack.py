import numpy as np 
import matplotlib.pylab as plt
from matplotlib import colors, cm
from astropy.io import fits, ascii
from astropy import units as u 
from astropy.cosmology import default_cosmology
import glob, os
from scipy import ndimage 

#plt.use_style('plot_style.txt')

lcdm = default_cosmology.get()

try:
	os.chdir('/data/gnarming/uchadaya/galreg/')
except FileNotFoundError:
	os.chdir('/Users/mila/Documents/Research/Postdoc/eFEDS-galaxies')

sdss = ascii.read('SDSS/sdss_output.csv',delimiter=',')
sdss['ssfr'] -= 9
ms = sdss['logMass']

def sort(sdss):
	sfr = sdss['ssfr']+sdss['logMass']
	sort = np.argsort(sfr)
	ind1 = int(len(sfr)/3)
	top = (sfr > sfr[sort[int(2*ind1)]])
	bottom = (sfr < sfr[sort[ind1]])
	return top, bottom


def stacked_image(files, lenx, leny, inputtype='file', size=15):
	img = np.zeros((lenx, leny))
	i = 0
	for file in files:
		if inputtype == 'file':
			f = fits.getdata(file)
		else:
			f = file
		x, y = f.shape
		# print(x,y)
		if y:
			starty = int((leny - y)/2)
			endy = starty+y
			startx = int((lenx - x)/2)
			endx = startx+x
			img[startx:endx,starty:endy] += f
	return ndimage.gaussian_filter(img, sigma=1, order=0)

def multiple_cbars(nrows, ncols, sharex, sharey):
	fig, ax = plt.subplots(nrows=nrows, ncols=ncols, sharex=sharex, sharey=sharey)
	plt.tight_layout(rect=[0,0,0.9,1],h_pad=1, w_pad=0.3)

	# im1 = plot(ax1) - make plot as needed
	# im2 = plot(ax2)
	# im3 = plot(ax3)

	#the following settings work for a plot with 3 rows
	cb_ax = fig.add_axes([0.85, 0.685, 0.04, 0.28]) # [left_edge, lower_edge, width, height]
	cbar1 = fig.colorbar(im1, cax=cb_ax) 
	cb_ax = fig.add_axes([0.85, 0.365, 0.04,0.28]) 
	cbar2 = fig.colorbar(im2, cax=cb_ax)
	cb_ax = fig.add_axes([0.85, 0.045, 0.04,0.28]) 
	cbar3 = fig.colorbar(im3, cax=cb_ax)
	return fig, ax, cbar1, cbar2, cbar3

def make_arrays(prefix, energy_range='_0520',suffix='',nbins=10, galsperfile=66,verbose=False): 
	files = glob.glob(prefix+'_profile%s*%s*fits' % (energy_range,suffix))
	files.sort()
	if verbose:
		print(files)
	cts = np.zeros((int(len(files)*galsperfile), nbins))
	exp = np.zeros((int(len(files)*galsperfile), nbins))
	area = np.zeros((int(len(files)*galsperfile), nbins))
	for i in range(len(files)):
			f = fits.getdata(files[i])
			counts = f['COUNTS']
			length = int(len(counts)/nbins)
			shape = (length, nbins)
			if verbose:
				print(shape)
			cts[i*length:(i+1)*length] = counts.reshape(shape)
			exp[i*length:(i+1)*length] = f['MEAN_SRC_EXP'].reshape(shape)
			area[i*length:(i+1)*length] = (f['AREA'] / (20**2)).reshape(shape) #to convert pixels to arcsec
			print(i)
	return cts, exp, area

#cross-matching is non-trivial because I had already split it into two groups
sdss1 = sdss[(ms > 10.2)*(ms < 10.7)][lvalid.astype(int)]
sdss2 = sdss[(ms > 10.7)*(ms < 11.2)][hvalid.astype(int)]

def format_plot(ax, yconv=3.1e39):
	ax.set_xlabel('R (kpc)')
	ax.set_xlim(4,100)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_ylim(3e33, 4e37)
	ax.set_ylabel(r'erg/s/kpc$^2$')
	for y in ax.get_yticklabels(): y.set_fontfamily('STIXGeneral')
	for x in ax.get_xticklabels(): x.set_fontfamily('STIXGeneral')
	plt.tight_layout()

from psf_contam import psf_curve, fit, xrb 
from psf_contam import area as area_psf
	
def profiles(cts, exp, area, sdss, yconv=3.1e39, simsfile = 'BenData/SBprof_disp.EAGLE_TNG.Ms10.20_10.70fixed.sSFR.docalc.20.list.SB_stats.out', 
	bckcols=5, ret=False, ax = None, xrb=False,cts_to_erg = 1.16e-12,color='tab:blue', ls='solid',beta=1.62e39, raw=False, log=False,verbose=False,bmin=None):
	 #cts/s to erg/s/cm**2

	#if we're only using T5 and T7, we only use photons at 0.8-2.0keV, but we still want Lx in 0.5-2.0, 
	#so instead the conversion factor is 2.26e-12. *7/2 since we now have only 2 out of 7 detectors 

	#similarly, if only using T12346, need to multiply cts_to_ergs by 7/5.

	z = sdss['z']
	dcorr = 4*np.pi*lcdm.luminosity_distance(z).to('cm')**2 #for luminosity distance
	dA = lcdm.angular_diameter_distance(z).to('kpc').value*u.arcsec.to('radian')
	
	area_kpc = np.zeros(area.shape)
	if log:
		xi = 10**np.linspace(0, np.log10(150),cts.shape[1]+1)
		x = (xi[1:]+xi[:-1])/2.
		xerr = (xi[1:] - xi[:-1])/2.
		if verbose:
			print(x)
	else:
		x = np.arange(0, cts.shape[1], 1)*10 + 5
		xerr = 5

	net_cts = np.zeros(cts.shape)
	if xrb:
		norm = Lxrb(sdss,beta=beta).mean()
		xp = np.arange(0,100,.1)
		yp = psf_curve(xp, *fit)/area_psf
		sbxrb = yp*norm
		fit_xrb = fit[0]*norm/area_psf
		lxrb = np.array([quad(psf_area, a=(x[i]-xerr[i]), b=(x[i]+xerr[i]),args=(fit_xrb, fit[1], fit[2]))[0] for i in range(len(x))]) 
			   #norm is erg/kpc**2 * kpc**2 = erg/s
		area_x = np.pi*((x+5)**2 - (x-5)**2) #in kpc^2
		sb_x = lxrb/area_x	
			  #this is in erg/s/kpc^2 
		xrb_cts_sec = lxrb / (cts_to_erg * np.mean(dcorr.value))
	else:
		xrb_cts_sec = np.zeros(x.shape)

	"""
	Plotting SB bkg in erg/s assumes that the background is coming from the same redshift as the galaxy
	MAYBE don't show background lines, because that is in counts
	"""

	rows = np.zeros(cts.shape)
	if bmin:
		src_cts = np.zeros(len(cts))
	for i in range(exp.shape[0]):
		area_kpc[i] = area[i]*(dA[i]**2)  #arcsec^2 to kpc^2
		bk_rate     = (cts[i,-bckcols:]/area[i,-bckcols:])
		src_cts[i]  = np.nansum(cts[i, bmin:-bckcols] - bk_rate*area[i,bmin:-bckcols])
		if xrb:
			rows[i] = cts[i] - xrb_cts_sec*exp[i]*area_kpc[i]/area_x
			if verbose:
				print("XRB subtracted from total profile, if requested")
			src_cts[i] -= np.nansum(xrb_cts_sec*exp[i,bmin:-bckcols]*area_kpc[i,bmin:-bckcols]/area_x)
		else:
			rows[i] = cts[i] 
		
		for j in range(exp.shape[1]):
			bkg = np.nanmean(rows[i,-bckcols:]/area[i,-bckcols:]) * area[i,j]
			if not exp[i,j]:
				exp[i,j] = np.nan
				cts[i,j] = np.nan
			net_cts[i,j] = rows[i,j] - bkg
	if verbose:
		print("bck computed from outer annuli")

	if raw:
		if bmin:
			return dcorr, src_cts
		else:
			return net_cts, area_kpc, dcorr

	else:
		sb = np.nansum(net_cts/area_kpc,axis=0) / np.nansum(exp, axis=0) 
		sb *= cts_to_erg * np.mean(dcorr.value)
		if verbose:
			print(sb)

		#background line
		sb_bkg = np.nanmean(rows[:,-bckcols:]/(area[:,-bckcols:]*exp[:,-bckcols:]))
		sb_bkg *= cts_to_erg * np.mean(dcorr.value)
		if verbose: 
			print("SB converted to erg/kpc**2")
		
		err = np.sqrt(np.nansum(cts, axis=0)) / np.nansum(net_cts, axis=0)
		errtot = np.sqrt(np.nansum(cts))/np.nansum((net_cts))

		if ax:
			ax.errorbar(x, sb, xerr=xerr, yerr=err*sb, linewidth=0, elinewidth=1, color=color)
			if simsfile:
				ax.hlines(sb_bkg, 0,200, color='k', alpha=0.25)
				ax.hlines(0.05*sb_bkg, 0,200, color='k', alpha=0.25, linestyle='dotted')
				sims = np.genfromtxt(simsfile)
				ax.fill_between(sims[:,0],sims[:,2]*yconv,sims[:,3]*yconv,color='tab:blue',alpha=0.5)
				ax.fill_between(sims[:,0],sims[:,5]*yconv,sims[:,6]*yconv,color='tab:red',alpha=0.5)
				ax.fill_between(sims[:,0],sims[:,8]*yconv,sims[:,9]*yconv,color='tab:purple',alpha=0.5)
				ax.fill_between(sims[:,0],sims[:,-2]*yconv,sims[:,-1]*yconv,color='tab:orange',alpha=0.5)
			if xrb:
				ax.plot(xp, sbxrb, c='k', label=r'XRB $\times$ PSF', linestyle=ls)
			format_plot(ax)
			return fig, ax 
		elif ret:
			return x, sb, errtot, sb_bkg	

#actually ya - I don't have to use up3 or down3! i can simply use the split function
from astropy.coordinates import SkyCoord
from astropy import units as u

def crossmatch(gal_region_file, sdss):
	#gal_region_file contains the dmcopy or dmextract regions
	#sdss is a table

	tra = [l.split('(')[1].split(',')[0] for l in gal_region_file]
	tdec = [l.split(',')[1] for l in gal_region_file]
	inds = []
	for i in range(len(tra)):
		radec = SkyCoord(tra[i]+' '+tdec[i], unit=(u.hourangle, u.deg))
		inds.append(np.argmin(abs(sdss['ra'] - radec.ra.to('deg').value) + abs(sdss['dec'] - radec.dec.to('deg').value)))
	return inds

def agn_in_gal(pts_cat, sdss, dmax=40):
	dmax *= u.arcsec.to('deg')
	tra = [l.split('(')[1].split(',')[0] for l in pts_cat]
	tdec = [l.split(',')[1] for l in pts_cat]
	radec = [SkyCoord(tra[i]+' '+tdec[i], unit=(u.hourangle, u.deg)) for i in range(len(tra))]
	radec = np.array([(radec[i].ra.to('deg').value, radec[i].dec.to('deg').value) for i in range(len(radec))])
	dist2 = np.zeros((len(sdss),len(pts_cat)))
	for i in range(len(sdss)):
		dist2[i] = [(abs(sdss['ra'][i] - radec[j,0])**2 + abs(sdss['dec'][i] - radec[j,1])**2) for j in range(len(radec))]
		print(i) #this is in degrees
	hasagn = [bool(sum(dist2[i] < dmax**2)) for i in range(len(sdss))]
	return dist2, hasagn


def Lx_stack(Mright, Mleft, cts, area, exp, sdss, bmin=1,bckcols=5, xrb=False, beta=1.62e39, cts_to_erg=1.624e-12):
	ms = sdss['logMass']
	Lx = np.zeros(len(Mright))
	binerr = np.zeros(len(Mright))

	for i in range(len(Mright)):
		# select galaxies in a given mass bin
		mask  = (ms > Mleft[i])*(ms < Mright[i])
		mask = [j for j in range(len(mask)) if mask[j]]
		scts = cts[mask]
		sexp = exp[mask]
		sz = sdss[mask]
		sarea = area[mask] #in arcsec
		dcorr = 4*np.pi*lcdm.luminosity_distance(sz['z']).to('cm')**2 #for luminosity distance

		bkg = np.nanmean(scts[:,-bckcols:]/(sarea[:,-bckcols:]*sexp[:,-bckcols:]))
		net_cts=np.nansum(scts[:,bmin:-bckcols] - bkg*sarea[:,bmin:-bckcols]*sexp[:,bmin:-bckcols])
		"Dec 28 - Added exposure weighting to BKG"

		Lx[i] = net_cts / np.nansum(sexp[:,bmin:-bckcols]) * cts_to_erg * dcorr.mean().value

		if xrb:
			norm = Lxrb(sz,beta=beta).mean()
			fit_xrb = fit[0]*norm/area_psf
			lxrb, _ = quad(psf_area, a=x[bmin], b=x[-bckcols],args=(fit_xrb, fit[1], fit[2]))
				   #norm is erg/kpc**2 * kpc**2 = erg/s
			print(Lx[i]/lxrb)	   
			xrb_cts_sec = lxrb *np.nanmean(sexp[:,bmin:-bckcols])/ (cts_to_erg * np.mean(dcorr.value))
			Lx[i] -= np.nansum(lxrb)
			net_cts -= lxrb*np.nansum(sexp[:,bmin:-bckcols])/(cts_to_erg * dcorr.mean().value)
		
		binerr[i] = np.sqrt(np.nansum(scts[:,bmin:-bckcols]))/net_cts
	return Lx, binerr

def Lx_split(Mright, Mleft, regionfile, cts, area, exp, sdss, bckcols=5, xrb=False, beta=1.62e39):	
		inds = crossmatch(regionfile, sdss)
		sdsssorted = sdss[inds]
		top, bottom = sort(sdsssorted)
		Lx_blue, binerr_blue = Lx_stack(Mright, Mleft, cts[top], area[top], exp[top], sdsssorted[top], condn = 'up3', xrb=xrb, beta=beta) 
		Lx_red, binerr_red = Lx_stack(Mright, Mleft, cts[bottom], area[bottom], exp[bottom], sdsssorted[bottom], condn = 'down3', xrb=xrb, beta=beta) 
		return Lx_blue, binerr_blue, Lx_red, binerr_red 
	
def Lx_basic_plot(Lx, binerr, Mstar, Mleft, sdss=sdss, ms=ms, xrb=True, fig=None, ax=None, beta=1.62e39, sims=False):	
	if not fig:
		fig, ax = plt.subplots()
	dM = Mstar - Mleft
	xerr = 10**Mstar - 10**Mleft

	Lxrb = np.zeros (len(Mstar))
	err_xrb = np.zeros ((2, len(Mstar)))
	plt.errorbar(10**Mstar, Lx, xerr=xerr, yerr=binerr*Lx, linewidth=0, elinewidth=1,color='tab:blue')
	plt.xscale('log')
	plt.yscale('log')
	plt.xticks(np.array([2e10,6e10,1e11]),[r'2$\times10^{10}$',r'6$\times10^{10}$',r'1$\times10^{11}$'],
				fontsize=14)
	plt.yticks(np.array([1e37,1e39,1e41,1e43]),fontsize=14)
	plt.xlim(10**10.15, 10**11.25)
	plt.ylim(5e36,1e42)
	plt.xlabel(r'$M_* / M_\odot$',fontsize=14,labelpad=0)
	plt.ylabel(r'$L_x$ (erg/s)',fontsize=14,labelpad=0)

	if sims:
		eagle0510 = np.vstack((np.genfromtxt('BenData/EAGLE_High-Mass.0.5_1.0keV.list.cat'), np.genfromtxt('BenData/EAGLE_Low-Mass.0.5_1.0keV.list.cat')))
		eagle1020 = np.vstack((np.genfromtxt('BenData/EAGLE_High-Mass.1.0_2.0keV.list.cat'), np.genfromtxt('BenData/EAGLE_Low-Mass.1.0_2.0keV.list.cat')))
		tng0510 = np.vstack((np.genfromtxt('BenData/TNG_High-Mass.0.5_1.0keV.list.cat'), np.genfromtxt('BenData/TNG_Low-Mass.0.5_1.0keV.list.cat')))
		tng1020 = np.vstack((np.genfromtxt('BenData/TNG_High-Mass.1.0_2.0keV.list.cat'), np.genfromtxt('BenData/TNG_Low-Mass.1.0_2.0keV.list.cat')))

		#colnum 3 = MStar, 4 = sSFR, -3 = LX_soft_ext
		plt.scatter(10**tng0510[:,3], (10**tng0510[:,-3] + 10**tng1020[:,-3]), alpha=0.1, color='tab:orange')
		plt.scatter(10**eagle0510[:,3], (10**eagle0510[:,-3]+10**eagle1020[:,-3]), alpha=0.1, color='darkgreen')

	for y in ax.get_yticklabels(): y.set_fontfamily('STIXGeneral')
	for x in ax.get_xticklabels(): x.set_fontfamily('STIXGeneral')
	
	return fig, ax 

def overplot_xrb(Mstar, fig, ax, sdss=sdss, ms=ms, beta=1.62e39):
	dM = (Mstar[1] - Mstar[0])/2.
	Lxrb = np.zeros(len(Mstar))
	err_xrb = np.zeros((2, len(Mstar)))
	for i in range(len(Mstar)):
		mask  = ((ms > (Mstar[i] - dM))*(ms < (Mstar[i] + dM)))
		mask = [j for j in range(len(mask)) if mask[j]]
		Lxrb[i] = Lxrb(sdss[mask], beta=beta).mean()
		Lxrb_low = Lxrb(sdss[mask], low=True,beta=beta).mean()
		Lxrb_high = Lxrb(sdss[mask], high=True,beta=beta).mean()
		err_xrb[0,i] = Lxrb[i] - Lxrb_low
		err_xrb[1,i] = Lxrb_high - Lxrb[i] 
	ax.errorbar(10**Mstar, Lxrb, xerr=xerr, yerr=err_xrb, linewidth=0, elinewidth=1,color='k')
	return Lxrb

def rebin(cts, net_cts, area_kpc, exp, dcorr, x, xmin, xmax, cts_to_erg = 1.16e-12*7/5):
	imin = np.argmin(abs(x-xmin))
	imax = np.argmin(abs(x-xmax))
	sub_cts = cts[:, imin:imax+1]
	sub_net = net_cts[:, imin:imax+1]
	sub_area = area_kpc[:, imin:imax+1]
	sub_exp = exp[:, imin:imax+1]
	sb = np.nansum(sub_net/sub_area)/np.nansum(sub_exp)
	sb *= cts_to_erg * np.mean(dcorr.value)
	err = np.sqrt(np.nansum(sub_cts))/np.nansum((sub_net))
	return sb, err

#combining mass bins for Lx-Mstar must be much easier, because I specify Mstar.

Why are the SF gals in the sims brighter? Is it because they have less massive AGN, and more AGN activity blows out the CGM? 
I.e. whatever quenches galaxies also blows out the CGM.

By eye estimate:
no XRB - low-M*  ~ 5e35*3.14*30**2 = 1.4e39 
	   - high-M* ~ 3.14e36*(281 + 131 + 163 + 1375+ 1282.5 + 1300) = 7.3e39 
no LXRB: low-M*  ~ 2e36*3.14*(400) = 2.5e39 - 2.5e38*pi = 1.7e39 
		 high-M* ~ 3.14e37*(56.25 + 24 + 37.5 + 86.25 + 180 + 142.5 + 156) = 2.1e40 - (2.53e39) = 1.85e40
raw:     low-M*  ~ 3.14e36*(6*56.25 + 4.5*(43.75) + 3.5*(125) + 3*(675) + 2125) = 1.65e40 - 5e38*3.14 = 1.5e40
		 high-M* ~ 3.14e37*(1019) = 3.2e40 - (2.65e39 + 1.37e39) = 2.8e40

Note that I have 3 mass bins in the current Lx-Mstar plot. So they will be a bit above/below these values
Though maybe I should try with two mass bins to test
