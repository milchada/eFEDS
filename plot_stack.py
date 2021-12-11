import numpy as np 
import matplotlib.pylab as plt
from matplotlib import colors, cm
from astropy.io import fits, ascii
from astropy import units as u 
from astropy.cosmology import default_cosmology
import glob, os
from scipy import ndimage 

plt.use_style('plot_style.txt')

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

def make_arrays(prefix, energy_range='_0520',nbins=10, galsperfile=66): 
	#'_0510', '' = '_0520', '1020' = '0520' - '0510'
	
	files = glob.glob(prefix+'_profile%s*fits' % energy_range) 
	files.sort()
	print(files)
	cts = np.zeros((int(len(files)*galsperfile), nbins))
	exp = np.zeros((int(len(files)*galsperfile), nbins))
	area = np.zeros((int(len(files)*galsperfile), nbins))
	ra = np.zeros(int(len(files)*galsperfile))
	dec = np.zeros(int(len(files)*galsperfile))
	for i in range(len(files)):
		# try:
			f = fits.getdata(files[i])
			counts = f['COUNTS']
			length = int(len(counts)/nbins)
			shape = (length, nbins)
			cts[i*length:(i+1)*length] = counts.reshape(shape)
			exp[i*length:(i+1)*length] = f['MEAN_SRC_EXP'].reshape(shape)
			area[i*length:(i+1)*length] = (f['AREA'] / (20**2)).reshape(shape) #to convert pixels to arcsec
			ra[i*length:(i+1)*length] = f['RA']
			dec[i*length:(i+1)*length] = f['DEC']
			print(i)
		# except (ValueError, KeyError):
		# 	print(i, ' oops!')
		# 	continue
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
	
def profiles(cts, exp, area, sdss, yconv=3.1e39, bck=None, simsfile = '../Ms10.20_10.70.out', 
	ret=True, bckcols=5, ax = None, xrb=False,cts_to_erg = 1.16e-12,color='tab:blue', ls='solid',beta=1.62e39):
	 #cts/s to erg/s/cm**2

	#if we're only using T5 and T7, we only use photons at 0.8-2.0keV, but we still want Lx in 0.5-2.0, 
	#so instead the conversion factor is 2.26e-12. *7/2 since we now have only 2 out of 7 detectors 

	#similarly, if only using T12346, need to multiply cts_to_ergs by 7/5.

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

	print("net counts computed")

	sb = np.nansum(net_cts/area_kpc,axis=0) / np.nansum(exp, axis=0) 
	#EVERYTHING SO FAR IS IN COUNTS
	sb *= cts_to_erg * np.mean(dcorr.value)
	
	#background line
	if bck:
		sb_bkg = np.nansum(bck['COUNTS']/(bck['AREA']*dA**2 / 400.))
		sb_bkg /= (np.nansum(bck['MEAN_SRC_EXP']))
	else:
		sb_bkg = np.nanmean(cts[:,-bckcols:]/(area[:,-bckcols:]*exp[:,-bckcols:]))
	sb_bkg *= cts_to_erg * np.mean(dcorr.value)

	print("SB converted to erg/kpc**2")
	
	#XRB contribution
	if xrb:
		norm = sb_xrb(sdss,beta=beta).mean()/area_psf[0]
		yp = psf_curve(x, *fit)
		sbxrb = yp*norm
		sb -= sbxrb
	"MODIFY THIS LINE SO THAT THE NET COUNT INCLUDES A CORRECTION FOR XRB"
	"THIS IS A FIRST PASS BUT ALMOST CERTAINLY INCOMPLETE"
		sbxrb /= cts_to_erg "times some area correction"
		for i in range(len(net_cts)):
			net_cts[i] -= sbxrb

	err = np.sqrt(np.nansum(cts, axis=0)) / np.nansum(net_cts, axis=0)
	errtot = np.sqrt(np.nansum(cts))/np.nansum((net_cts))

	if ax:
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
			ax.plot(x, y*norm, c='k', label=r'XRB $\times$ PSF', linestyle=ls)
		format_plot(ax)
		return fig, ax 
	elif ret:
		return x, sb, errtot, sb_bkg	

def Lx_stack(Mstar, cts, sdss, area, condn=None, exp=None, bmin=1,bckcols=5, xrb=False, beta=1.62e39):
	dM = (Mstar[1] - Mstar[0])/2.
	ms = sdss['logMass']
	Lx = np.zeros(len(Mstar))
	binerr = np.zeros(len(Mstar))

	for i in range(len(Mstar)):
		# select galaxies in a given mass bin
		mask  = (ms > (Mstar[i] - dM))*(ms < (Mstar[i] + dM))
		mask = [j for j in range(len(mask)) if mask[j]]

		if sum(mask):			
			scts = cts[mask]
			sexp = exp[mask]
			sz = sdss[mask]
			sarea = area[mask]
			print(len(scts))

			x, sb, binerr[i], sb_bkg = profiles(scts, sexp, sarea,sz, bckcols=bckcols, simsfile=None, xrb=xrb, beta=beta)
			# print(sb)
			x1 = x+5
			x0 = x-5
			ar = np.pi*(x1**2 - x0**2)
			
			ypos = (sb*ar)[:-bckcols] 
			Lx[i] = sum(ypos[ypos>0]) 
	return Lx, binerr

#actually ya - I don't have to use up3 or down3! i can simply use the split
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

def Lx_obs(regionfile, cts, area, exp, sdss, split = False, dM = 0.1, bckcols=5, xrb=False, beta=1.62e39):	
	Mstar = np.arange(10.2+dM,11.2,2*dM)
	if split:
		inds = crossmatch(regionfile, sdss)
		sdsssorted = sdss[inds]
		top, bottom = sort(sdsssorted)
		Lx_blue, binerr_blue = Lx_stack(Mstar, cts[top], sdsssorted[top], area[top], condn = 'up3', exp=exp[top], xrb=xrb, beta=beta) 
		Lx_red, binerr_red = Lx_stack(Mstar, cts[bottom], sdsssorted[bottom], area[bottom], condn = 'down3', exp=exp[bottom], xrb=xrb, beta=beta) 
		return Lx_blue, binerr_blue, Lx_red, binerr_red 
	else:
		Lx, binerr = Lx_stack(Mstar, cts, sdss, area, condn=None, exp=exp, xrb=xrb, beta=beta) 
		return Lx, binerr
	
def plot(Lx, binerr, Mstar, sdss=sdss, ms=ms, xrb=True, fig=None, ax=None, beta=1.62e39):	
	if not fig:
		fig, ax = plt.subplots()
	dM = (Mstar[1] - Mstar[0])/2.
	Mleft = Mstar - dM   
	xerr = 10**Mstar - 10**Mleft

	Lxrb = np.zeros (len(Mstar))
	err_xrb = np.zeros ((2, len(Mstar)))
	if xrb:
		for i in range(len(Mstar)):
			mask  = ((ms > (Mstar[i] - dM))*(ms < (Mstar[i] + dM)))
			mask = [j for j in range(len(mask)) if mask[j]]
			Lxrb[i] = sb_xrb(sdss[mask], beta=beta).mean()
			Lxrb_low = sb_xrb(sdss[mask], low=True,beta=beta).mean()
			Lxrb_high = sb_xrb(sdss[mask], high=True,beta=beta).mean()
			err_xrb[0,i] = Lxrb[i] - Lxrb_low
			err_xrb[1,i] = Lxrb_high - Lxrb[i] 
		plt.errorbar(10**Mstar, Lxrb, xerr=xerr, yerr=err_xrb, linewidth=0, elinewidth=1,color='k')
	plt.errorbar(10**Mstar, Lx, xerr=xerr, yerr=binerr*Lx, linewidth=0, elinewidth=1,color='tab:blue')
	plt.xscale('log')
	plt.yscale('log')

	return fig, ax 

def Lx_sims(fig, ax, split=False):
	Mstar = np.array([10.3, 10.5, 10.7, 10.9, 11.1])
	eagle = np.vstack((np.genfromtxt('../EAGLE_High-Mass.cat'), np.genfromtxt('../EAGLE_Low-Mass.cat')))
	tng = np.vstack((np.genfromtxt('../TNG_High-Mass.cat'), np.genfromtxt('../TNG_Low-Mass.cat')))
	#colnum 3 = MStar, 4 = sSFR, -3 = LX_soft_ext
	plt.scatter(10**tng[:,3], 10**tng[:,-3], alpha=0.1, color='tab:orange')
	plt.scatter(10**eagle[:,3], 10**eagle[:,-3], alpha=0.1, color='tab:red')


"OK how the fuck would I combine just the outer few bins in the SB plot"
if low:
	x0 = [0, 10, 20, 30, 40]
	x1 = [10,20, 30, 40, 100]
if high:
	x0 = [0, 10, 20, 30, 40, 50, 60]
	x0 = [10,20, 30, 40, 50, 60, 100]
x = (x1 - x0)/2.