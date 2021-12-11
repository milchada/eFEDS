import numpy as np
import matplotlib.pylab as plt
import glob, os
from astropy.io import fits, ascii

plt.style.use('seaborn-deep')
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Ubuntu'
plt.rcParams['font.monospace'] = 'Ubuntu Mono'
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['figure.titlesize'] = 14

lowbck  = fits.getdata('lowmass_bck_0520.fits')
lvalid  = np.load('lvalid.npy').astype(int)
lowbck  = lowbck[lvalid]

highbck = fits.getdata('highmass_bck_0520.fits')
# hvalid = [j for j in range(len(hcts)) if np.nansum(hexp[j])]
hvalid  = np.load('hvalid.npy').astype(int)
highbck = highbck[hvalid] 

cts    = np.load('cts0520.npy') 
exp    = np.load('exp0520.npy')
area   = np.load('area0520.npy')

try:
	sdss = ascii.read('../sdss_output.csv')
except FileNotFoundError:
	sdss = ascii.read('SDSS/sdss_output.csv')
sdss['ssfr'] -= 9
ms = sdss['logMass']

#cross-matching is non-trivial because I had already split it into two groups
sdss1 = sdss[(ms > 10.2)*(ms < 10.7)]#[lvalid.astype(int)]
sdss2 = sdss[(ms > 10.7)*(ms < 11.2)]#[hvalid.astype(int)]

def fig1_old():
	#obs not split by SFR

	fig, ax = plt.subplots(ncols=2, sharex=True, sharey=True)
	ax1, ax2 = ax.flatten()
	profiles(cts[:1233], exp[:1233], area[:1233], sdss1, yconv = 3.1e39, simsfile = '../Ms10.20_10.70.out', 
	ret_im=False, ret = False, bckcols=5, ax = ax1, xrb = True)
	profiles(cts[1233:], exp[1233:], area[1233:], sdss2, yconv = 3.1e39, simsfile = '../Ms10.70_11.20.out', 
	ret_im=False, ret = False, bckcols=5, ax = ax2, xrb = True)
	return fig, ax

def sort(sdss):
	sfr = sdss['ssfr']+sdss['logMass']
	sort = np.argsort(sfr)
	ind = int(len(sfr)/3)
	top = (sfr > sfr[sort[int(2*ind)]])
	bottom = (sfr < sfr[sort[ind]])
	return top, bottom

def fig1():
	fig, ax = plt.subplots(ncols=2, sharex=True, sharey=True)
	ax1, ax2 = ax.flatten()
	#split each mass bin into SF and quenched
	top, bottom = sort(sdss1)
	profiles(cts[:1233][top], exp[:1233][top], area[:1233][top], sdss1[top], yconv = 3.1e39, simsfile = 'BenData/Ms10.20_10.70.out', 
	ret_im=True, ret = False, bckcols=5, ax = ax1, xrb = True,color='teal')
	profiles(cts[:1233][bottom], exp[:1233][bottom], area[:1233][bottom], sdss1[bottom], yconv = 3.1e39, simsfile=None, 
	ret_im=True, ret = False, bckcols=5, ax = ax1, xrb = True,color='tab:red')
	top, bottom = sort(sdss2)
	profiles(cts[1233:][top], exp[1233:][top], area[1233:][top], sdss2[top], yconv = 3.1e39, simsfile = 'BenData/Ms10.70_11.20.out', 
	ret_im=True, ret = False, bckcols=5, ax = ax2, xrb = True,color='teal')
	profiles(cts[1233:][bottom], exp[1233:][bottom], area[1233:][bottom], sdss2[bottom], yconv = 3.1e39, simsfile=None,
	ret_im=True, ret = False, bckcols=5, ax = ax2, xrb = True,color='tab:red')
	return fig, ax


def fig2(dM = 0.1, xrb=True, sims=True, split=False):	
	
	Lx_blue, binerr_blue, Lx_red, binerr_red = Lx_obs(cts, area, exp, sdss, split = 2, dM = dM, bckcols=5)	
	Lx, binerr = Lx_obs(cts, area, exp, sdss, split = False, dM = dM, bckcols=5)	

	fig, ax = plt.subplots()
	Mstar = np.arange(10.2+dM,11.2,2*dM)
	Mleft = Mstar - dM   
	xerr = 10**Mstar - 10**Mleft

	Lxrb = np.zeros (len(Mstar))
	err_xrb = np.zeros ((2, len(Mstar)))
	if xrb:
		for i in range(len(Mstar)):
			mask  = ((ms > (Mstar[i] - dM))*(ms < (Mstar[i] + dM)))
			mask = [j for j in range(len(mask)) if mask[j]]
			Lxrb[i] = sb_xrb(sdss[mask]).mean()
			Lxrb_low = sb_xrb(sdss[mask], low=True).mean()
			Lxrb_high = sb_xrb(sdss[mask], high=True).mean()
			err_xrb[0,i] = Lxrb[i] - Lxrb_low
			err_xrb[1,i] = Lxrb_high - Lxrb[i] 
		plt.errorbar(10**Mstar, Lxrb, xerr=xerr, yerr=err_xrb, linewidth=0, elinewidth=1,color='k')
	plt.errorbar(10**Mstar, Lx_blue, xerr=xerr, yerr=binerr_blue*Lx_blue, linewidth=0, elinewidth=1,color='tab:blue')
	plt.errorbar(10**Mstar, Lx_red, xerr=xerr, yerr=binerr_red*Lx_red, linewidth=0, elinewidth=1,color='tab:red')
	plt.errorbar(10**Mstar, Lx, xerr=xerr, yerr=binerr*Lx, linewidth=0, elinewidth=1,color='tab:blue')
	plt.xscale('log')
	plt.yscale('log')

	if sims:
		eagle = np.vstack((np.genfromtxt('BenData/EAGLE_High-Mass.cat'), np.genfromtxt('BenData/EAGLE_Low-Mass.cat')))
		 tng0510 = np.vstack((np.genfromtxt('BenData/TNG_High-Mass.0.5_1.0keV.list.cat'), np.genfromtxt('BenData/TNG_Low-Mass.0.5_1.0keV.list.cat')))
		 tng0510 = np.vstack((np.genfromtxt('BenData/TNG_High-Mass.1.0_2.0keV.list.cat'), np.genfromtxt('BenData/TNG_Low-Mass.1.0_2.0keV.list.cat')))
		 tng = tng0510
		 tng[:,-3] = np.log10(10**tng0510[:,-3]+10**tng1020[:,-3])
		#colnum 3 = MStar, 4 = sSFR, -3 = LX_soft_ext
		plt.scatter(10**tng[:,3], 10**tng[:,-3], alpha=0.1, color='indigo')
		plt.scatter(10**eagle[:,3], 10**eagle[:,-3], alpha=0.1, color='silver')

	plt.xlabel(r'log(M$_*$/M$_\odot$)')
	plt.ylabel(r'L$_{X, 0.5-2.0 keV}$ (> 10 kpc)')
	plt.tight_layout()
	return fig, ax 

def fig3():
	#galaxy sample properties
	fig, ax = plt.subplots(ncols=3)
	ax1, ax2, ax3 = ax.flatten()
	hist1, bins1 = np.histogram(sdss['logMass'], range=(10.2,11.2), bins=20)
	hist2, bins2 = np.histogram((sdss['logMass']+sdss['ssfr']), range=(-18,3), bins=21)
	hist3, bins3 = np.histogram(sdss['z'], range=(0.01,0.1), bins=11)
	dx1 = (bins1[1] - bins1[0])/2.
	dx2 = (bins2[1] - bins2[0])/2.
	dx3 = (bins3[1] - bins3[0])/2.
	ax1.step(bins1[1:]-dx1, hist1, color='tab:blue')
	ax2.step(bins2[1:]-dx2, hist2, color='tab:blue')
	ax3.step(bins3[1:]-dx3, hist3, color='tab:blue')
