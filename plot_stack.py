from psf_contam import psf_curve, fit, xrb 
from psf_contam import area as area_psf
import numpy as np 
import matplotlib.pylab as plt
from matplotlib import colors, cm
from astropy.io import fits, ascii
from astropy import units as u 
from astropy.cosmology import default_cosmology
import glob, os
from scipy import ndimage 
from astropy.coordinates import SkyCoord
from astropy import units as u

lcdm = default_cosmology.get()

try:
	os.chdir('/data/gnarming/uchadaya/galreg/')
	sdss = ascii.read('sdss_output.csv',delimiter=',')
except FileNotFoundError:
	os.chdir('/Users/mila/Documents/Research/Postdoc/eFEDS-galaxies')
	sdss = ascii.read('SDSS/sdss_output.csv',delimiter=',')

def sort(sdss):
    q = np.argsort(sdss['ssfr'])[:int(len(sdss)/3.)]
    sf = np.argsort(sdss['ssfr'])[int(2*len(sdss)/3.):]
    try: 
        return sf.values, q.values
    except AttributeError:
        return sf, q

def annulus(array, rmin, rmax):
    x = np.arange(array.shape[0])
    y = np.arange(array.shape[1])
    x = x - x.mean()
    y = y - y.mean()
    X, Y = np.meshgrid(x,y)
    r2 = X**2 + Y**2
    return np.ma.masked_inside(r2, rmin**2, rmax**2).mask

def stacked_image(files, lenx, leny, inputtype='file'):
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

def stack_sf_q(dir, mean_exp, prefix='low', rmin=50, rmax=100, bckmin=25, bckmax=37.5):
	os.chdir(dir)
	gals = ascii.read('gal_props.csv', delimiter=',',header_start=0,)
	quiescent = gals['ind'][np.argsort(gals['ssfr'])[:int(len(gals)/3.)]]
	active = gals['ind'][np.argsort(gals['ssfr'])[int(2*len(gals)/3.):]]
	qfiles = [prefix+'mass_cutout_%d.fits' % ind for ind in quiescent]
	sffiles = [prefix+'mass_cutout_%d.fits' % ind for ind in active]
	im_q = stacked_image(qfiles, 151, 151)[rmin:rmax, rmin:rmax]
	im_sf = stacked_image(sffiles, 151, 151)[rmin:rmax, rmin:rmax]
	im_q /= (mean_exp*len(qfiles))
	im_sf /= (mean_exp*len(qfiles))
	bck_ind = annulus(im_q, bckmin, bckmax)
	q_bck = np.nanmedian(im_q[bck_ind])
	sf_bck = np.nanmedian(im_sf[bck_ind])
	return im_q - q_bck, im_sf - sf_bck

def plot_sf_q(rmin=50,rmax=100,bckmin=25, bckmax=37.5,c=62.5):
	fig, ax = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
	ax1, ax2, ax3, ax4 = ax.flatten()
	lexp = np.nanmean(np.load('/data/gnarming/uchadaya/galreg/lexp_0520.npy')/5., axis=1)
	lexp = np.nanmean(lexp[lexp > 1])
	hexp = np.nanmean(np.load('/data/gnarming/uchadaya/galreg/hexp_0520.npy')/5., axis=1)
	hexp = np.nanmean(hexp[hexp > 1])
	imlo_q, imlo_sf = stack_sf_q('/data/gnarming/uchadaya/galreg/lowmass/cutouts/0520/', lexp,rmin=rmin,rmax=rmax,bckmin=bckmin, bckmax=bckmax)
	imhi_q, imhi_sf = stack_sf_q('/data/gnarming/uchadaya/galreg/highmass/cutouts/0520/', hexp, prefix='high',rmin=rmin,rmax=rmax,bckmin=bckmin, bckmax=bckmax)
	zmax = max(imlo_q.max(), imlo_sf.max(), imhi_q.max(), imhi_sf.max())
	ax1.imshow(imlo_sf, origin='lower', cmap=cm.RdBu_r, norm=colors.Normalize(-zmax, zmax))
	ax2.imshow(imhi_sf, origin='lower', cmap=cm.RdBu_r, norm=colors.Normalize(-zmax, zmax))
	ax3.imshow(imlo_q, origin='lower', cmap=cm.RdBu_r, norm=colors.Normalize(-zmax, zmax))
	ax4.imshow(imhi_q, origin='lower', cmap=cm.RdBu_r, norm=colors.Normalize(-zmax, zmax))
	for a in ax.flatten():
		src_circle = plt.Circle((c, c), bckmin, fill=False, color='k', linestyle='dotted')
		bck_circle = plt.Circle((c, c), bckmax, fill=False, color='k', linestyle='dotted')
		a.add_patch(src_circle)
		a.add_patch(bck_circle)
    #for paper: fig, ax = plot_sf_q(13,138,25,62.5)
	return fig, ax 

def cts_area_ratio(tab, prof=False):
	if prof:
		cts = tab['cts-0-5.5'].sum()
		bck_cts = tab['bck_cts'].sum()
		area = tab['bck_area']/tab['area-0-5.5']
	else:
		cts = tab['cts'].sum()
		bck_cts = tab['bck-cts'].sum()
		area = tab['bck-area']/tab['bck-cts']
	return cts, bck_cts, np.nanmean(area[area>0])

def hardness_ratio(soft, hard, prof=False):
	softsrc, softbck, softarea = cts_area_ratio(soft, prof=prof)
	allsrc, allbck, allarea = cts_area_ratio(hard, prof=prof) 
	print(softsrc, softbck, softarea, (allsrc-softsrc), (allbck-softbck), allarea)
	return softsrc, softbck, softarea, (allsrc-softsrc), (allbck-softbck), allarea

def plot_hardness_ratios(xerr=0.25):
	def med_min_max(file):
		with open(file) as f:
			a = f.readlines()
		for line in a:
			if '(H-S)/(H+S)' in line:
				_,_,_,med, min, max,_ = line.split('	')
				return float(med), float(med)-float(min), float(max)-float(med)

	x = [(10.7+11.2)/2., (10.2+10.7)/2.]
	sf_core = np.zeros((2,3))
	sf_all = np.zeros((2,3))
	q_core = np.zeros((2,3))
	q_all = np.zeros((2,3))
	mean_core = np.zeros((2,3))
	mean_all = np.zeros((2,3))

	fig, ax = plt.subplots(ncols=2, sharex=True, sharey=True)
	
	q_core[1] = med_min_max('highcore_q.txt')
	sf_core[1] = med_min_max('highcore_sf.txt')
	mean_core[1] = med_min_max('highmass-0-10.txt')  

	mean_all[1] = med_min_max('highmass-10-100.txt')
	sf_all[1] = med_min_max('highsrc_sf.txt')
	q_all[1] = med_min_max('highsrc_q.txt')
	
	
	mean_core[0] = med_min_max('lowmass-0-10.txt')  
	sf_core[0] = med_min_max('lowcore_sf.txt')
	q_core[0] = med_min_max('lowcore_q.txt')

	mean_all[0] = med_min_max('lowmass-10-100.txt')  
	sf_all[0] = med_min_max('lowsrc_sf.txt')
	q_all[0] = med_min_max('lowsrc_q.txt')

	ax[0].errorbar(x, q_core[:,0], xerr=xerr, yerr=q_core[:,1:], color='tab:red', linewidth=0, elinewidth=2)
	ax[0].errorbar(x, sf_core[:,0], xerr=xerr, yerr=sf_core[:,1:], color='tab:blue', linewidth=0, elinewidth=2)
	ax[0].errorbar(x, mean_core[:,0], xerr=xerr, yerr=mean_core[:,1:], color='k', linewidth=0, elinewidth=2)
	
	ax[1].errorbar(x, q_all[:,0], xerr=xerr, yerr=q_all[:,1:], color='tab:red', linewidth=0, elinewidth=2)
	ax[1].errorbar(x, mean_all[:,0], xerr=xerr, yerr=mean_all[:,1:], color='k', linewidth=0, elinewidth=2)
	ax[1].errorbar(x, sf_all[:,0], xerr=xerr, yerr=sf_all[:,1:], color='tab:blue', linewidth=0, elinewidth=2)
	return fig, ax

	
def axisticks(ax, tfile, xmin, xmax, x=True, y=True, tickmax=400, ntix=5, factor=1, fontsize=14):
	xtix = np.linspace(-tickmax, tickmax, ntix)
	dx = fits.getheader(tfile)['CDELT1']*factor
	xpos = xtix/dx + (xmin+xmax)/2.
	if x:
		ax.set_xticks(xpos)
		ax.set_xticklabels(['%d' % x for x in xtix], fontsize=fontsize)
	if y:
		ax.set_yticks(xpos)
		ax.set_yticklabels(['%d' % x for x in xtix], fontsize=fontsize)

def make_arrays(prefix, energy_range='_0520',suffix='',nbins=10, galsperfile=66,verbose=False): 
	files = glob.glob(prefix+'profile%s*%s*fits' % (energy_range,suffix))
	files.sort()
	if verbose:
		print(files)
	cts = np.zeros((int(len(files)*galsperfile), nbins))
	exp = np.zeros((int(len(files)*galsperfile), nbins))
	area = np.zeros((int(len(files)*galsperfile), nbins))
	radec = np.zeros((int(len(files)*galsperfile), 2))
	for i in range(len(files)):
			f = fits.getdata(files[i])
			counts = f['COUNTS']
			length = int(len(counts)/nbins)
			shape = (length, nbins)
			if verbose:
				print(shape)
			cts[i*length:(i+1)*length] = counts.reshape(shape)
			try:
				exp[i*length:(i+1)*length] = f['MEAN_SRC_EXP'].reshape(shape)
			except KeyError:
				print('')
			area[i*length:(i+1)*length] = (f['AREA'] / (20**2)).reshape(shape) #to convert pixels to arcsec
			# print(f['RA'])
			ngal = len(f['RA'][::nbins]) #to convert pixels to arcsec
			radec[i*galsperfile:(i*galsperfile)+ngal,0] = f['RA'][::nbins] #to convert pixels to arcsec
			radec[i*galsperfile:(i*galsperfile)+ngal,1] = f['DEC'][::nbins] #to convert pixels to arcsec
			print(i)
	return cts, exp, area, radec

import pandas as pd
low = pd.read_csv('lowmass_profiles.csv')
high = pd.read_csv('highmass_profiles.csv')

def format_plot(ax, yconv=3.1e39,y=True):
	ax.set_xlabel('R (kpc)', fontsize=14)
	ax.set_xlim(5.5,100)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_ylim(3e33, 4e37)
	if y:
		ax.set_ylabel(r'$\Sigma_X$ (erg/s/kpc$^2$)',fontsize=14)
	for y in ax.get_yticklabels(): y.set_fontfamily('STIXGeneral')
	for x in ax.get_xticklabels(): x.set_fontfamily('STIXGeneral')
	plt.tight_layout()

def plot_obs(table, ax, color):
	tva = table[abs(table[:,-1]) < 1]
	ax.errorbar((tva[:,0]+ tva[:,1])/2., tva[:,2],xerr=(tva[:,0]- tva[:,1])/2., 
		yerr=tva[:,2]*tva[:,3], linewidth=0, elinewidth=1, color=color)
	ax.scatter((tva[:,0]+tva[:,1])/2., tva[:,2],color=color)

	tinva = table[abs(table[:,-1]) > 1]
	ax.scatter((tinva[:,0]+tinva[:,1])/2., tinva[:,2],color=color, marker='v')

def xrb_psf(sdss, beta=1.62e39, fit=fit, area_psf=area_psf):
		norm = Lxrb(sdss,beta=beta).mean()
		xp = np.arange(0,100,.1)
		yp = psf_curve(xp, *fit)/area_psf
		return xp, yp,norm

def profiles(profile_csv, bckcols=5, ret=False, xrb=False, cts_to_erg = 1.624e-12, beta=1.62e39, log=False,bmin=None):
	z = profile_csv['z']
	dcorr = 4*np.pi*lcdm.luminosity_distance(z).to('cm')**2 #for luminosity distance
	dA = lcdm.angular_diameter_distance(z).to('kpc').value*u.arcsec.to('radian')

	cts = [profile_csv[col] for col in profile_csv.columns if 'cts' in col]
	cts = np.array([l.values for l in cts]).T

	if log:
		xi = 10**np.linspace(0, np.log10(150),cts.shape[1]+1)
		x = (xi[1:]+xi[:-1])/2.
		xerr = (xi[1:] - xi[:-1])/2.
		area_x = np.pi*(xi[1:]**2 - xi[:-1]**2)
	else:
		x = np.arange(0, cts.shape[1], 1)*10 + 5
		xerr = 5
		area_x = np.pi*((x+5)**2 - (x-5)**2) #in kpc^2
	if xrb:
		xp, yp, norm = xrb_psf(profile_csv, beta=beta)
		sbxrb = yp*norm
		fit_xrb = fit[0]*norm/area_psf
		lxrb = np.array([quad(psf_area, a=(x[i]-xerr[i]), b=(x[i]+xerr[i]),args=(fit_xrb, fit[1], fit[2]))[0] for i in range(len(x))]) 
			   #norm is erg/kpc**2 * kpc**2 = erg/s
		sb_x = lxrb/area_x	
			  #this is in erg/s/kpc^2 
		xrb_cts_sec = lxrb / (cts_to_erg * np.mean(dcorr.value))
	else:
		xrb_cts_sec = np.zeros(x.shape)

	net_cts = np.zeros(cts.shape)
	sb_bkg = np.nanmedian(profile_csv['cts_14']/(profile_csv['area_14']* profile_csv['exp_14']))
	for j in range(cts.shape[1]):
		exp = profile_csv['exp_%d' % j]
		counts = profile_csv['cts_%d' % j]
		area = profile_csv['area_%d' % j]
		area_kpc = area*(dA**2)
		row = counts - xrb_cts_sec[j]*exp*area_kpc/area_x[j]
		net_cts[:,j] = row - sb_bkg*exp*area
        
	sb_bkg *= cts_to_erg * np.mean(dcorr.value)
	exp = [profile_csv[col] for col in profile_csv.columns if 'exp' in col]
	exp = np.array([l.values for l in exp]).T
	area_kpc = [profile_csv[col]*(dA**2) for col in profile_csv.columns if 'area' in col]
	area_kpc = np.array([l.values for l in area_kpc]).T

	return cts, area_kpc, exp, net_cts, dcorr, sb_bkg

def plot_sim_profiles(cts, net_cts, area_kpc, dcorr, exp, ax, simsfile = 'BenData/SBprof_disp.EAGLE_TNG.Ms10.20_10.70fixed.sSFR.docalc.20.list.SB_stats.out', 
	cts_to_erg = 1.624e-12,color='tab:blue', ls='solid',bckcols=1, xrb=False, beta=1.62e39):

	sb = np.nansum(net_cts/area_kpc,axis=0) / np.nansum(exp, axis=0) 
	sb *= cts_to_erg * np.mean(dcorr.value)

	#background line
	sb_bkg = np.nanmedian(rows[:,-bckcols:]/(area[:,-bckcols:]*exp[:,-bckcols:]))
	sb_bkg *= cts_to_erg * np.mean(dcorr.value)
	
	err = np.sqrt(np.nansum(cts, axis=0)) / np.nansum(net_cts, axis=0)
	errtot = np.sqrt(np.nansum(cts))/np.nansum((net_cts))

	ax.hlines(sb_bkg, 0,200, color='k', alpha=0.25)
	ax.hlines(0.05*sb_bkg, 0,200, color='k', alpha=0.25, linestyle='dotted')
	
	if simsfile:
		sims = np.genfromtxt(simsfile)
		yconv=3.1e39
		ax.fill_between(sims[:,0],sims[:,2]*yconv,sims[:,3]*yconv,color='tab:blue',alpha=0.5)
		ax.fill_between(sims[:,0],sims[:,5]*yconv,sims[:,6]*yconv,color='tab:red',alpha=0.5)
		ax.fill_between(sims[:,0],sims[:,8]*yconv,sims[:,9]*yconv,color='tab:purple',alpha=0.5)
		ax.fill_between(sims[:,0],sims[:,-2]*yconv,sims[:,-1]*yconv,color='tab:orange',alpha=0.5)
		format_plot(ax, y=y)

	if xrb:
		xp, yp, norm = xrb_psf(sdss, beta=beta)
		sbxrb = yp*norm
		ax.plot(xp, sbxrb, c='k', label=r'XRB $\times$ PSF', linestyle=ls)

def crossmatch(gal_region_file, sdss, tab=False, arr=False):
	#gal_region_file contains the dmcopy or dmextract regions
	#sdss is a table
	inds = []
	if tab:
		for i in range(len(gal_region_file)):
			inds.append(np.argmin(abs(sdss['ra'] - gal_region_file['ra'].iloc[i]) + abs(sdss['dec'] - gal_region_file['dec'].iloc[i])))
	elif arr:
		for i in range(len(gal_region_file)):
			inds.append(np.argmin(abs(sdss['ra'] - gal_region_file[i,0]) + abs(sdss['dec'] - gal_region_file[i,1])))
	else:
		with open(gal_region_file) as f:
			a = f.readlines()
		c = a[::15]
		print(len(a), len(b), len(c))
		tra = [l.split('(')[1].split(',')[0] for l in c]
		tdec = [l.split(',')[1] for l in c]
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

def Lx_stack(Mright, Mleft, src, rmin=10,rmax=100, xrb=False, beta=1.62e39, cts_to_erg=1.624e-12, 
             fit=fit, area_psf = area_psf):
	ms = src['logMass']
	Lx = np.zeros(len(Mright))
	binerr = np.zeros(len(Mright))

	for i in range(len(Mright)):
		# select galaxies in a given mass bin
		mask  = (ms > Mleft[i])*(ms < Mright[i])
		mask = [j for j in range(len(mask)) if mask[j]]
		table = src.iloc[mask]
		dcorr = 4*np.pi*lcdm.luminosity_distance(table['z']).to('cm')**2 #for luminosity distance

		bkg_rt = np.nanmedian(table['bck_cts']/(table['bck_area']*table['bck_exp']))
		net_cts=table['cts'] - bkg_rt*table['area']*table['exp']
		Lx[i] = np.nansum(net_cts * dcorr.value) / np.nansum(table['exp']) * cts_to_erg
# 		print(np.sqrt(np.nansum(table['cts']))/np.nansum(net_cts))
		if xrb:
			norm = Lxrb(table,beta=beta).mean()
			fit_xrb = fit[0]*norm/area_psf
			lxrb, _ = quad(psf_area, a=rmin, b=rmax,args=(fit_xrb, fit[1], fit[2]))
			xrb_cts_sec = lxrb *np.nanmean(table['exp'])/(cts_to_erg * np.mean(dcorr.value))
			Lx[i] -= np.nansum(lxrb)
			net_cts -= lxrb*np.nanmean(table['exp'])/(cts_to_erg * dcorr.mean().value)
		binerr[i] = np.sqrt(np.nansum(table['cts']))/np.nansum(net_cts)
	return Lx, binerr

def Lx_split(Mright, Mleft, regionfile, cts, area, exp, sdss, bckcols=5, xrb=False, beta=1.62e39):	
		inds = crossmatch(regionfile, sdss)
		sdsssorted = sdss[inds]
		top, bottom = sort(sdsssorted)
		Lx_blue, binerr_blue = Lx_stack(Mright, Mleft, cts[top], area[top], exp[top], sdsssorted[top], condn = 'up3', xrb=xrb, beta=beta) 
		Lx_red, binerr_red = Lx_stack(Mright, Mleft, cts[bottom], area[bottom], exp[bottom], sdsssorted[bottom], condn = 'down3', xrb=xrb, beta=beta) 
		return Lx_blue, binerr_blue, Lx_red, binerr_red 
	
def Lx_basic_plot(fig=None, ax=None):	
	if not fig:
		fig, ax = plt.subplots()
	plt.xscale('log')
	plt.yscale('log')
	plt.xticks(np.array([2e10,6e10,1e11]),[r'2$\times10^{10}$',r'6$\times10^{10}$',r'1$\times10^{11}$'],
				fontsize=14)
	plt.yticks(np.array([1e37,1e39,1e41,1e43]),fontsize=14)
	plt.xlim(10**10.15, 10**11.25)
	plt.ylim(5e36,1e42)
	plt.xlabel(r'$M_* / M_\odot$',fontsize=14,labelpad=0)
	plt.ylabel(r'$L_x$ (erg/s)',fontsize=14,labelpad=0)

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
