from astropy.io import fits
from astropy.coordinates import SkyCoord
import glob, os
import numpy as np

basedir = '/data/gnarming/uchadaya/galreg/lowmass/cutouts/0520'
os.chdir(basedir)
files = glob.glob('*cutout*fits'); files.sort()
def ra_dec(files):
	ras = np.zeros(len(files))
	decs = np.zeros(len(files))
	for i in range(len(files)):
		file = files[i]
		h = fits.getheader(file)['HISTORY']
		for l in range(len(h)):
			if 'box' in h[l]:
				print(h[l])
				ra = h[l].split('(')[-1].split('A')[0]
				ra += h[l+1].split(' :')[1].split(',')[0]
				dec = h[l+1].split(' :')[1].split(',')[1]
				break
		c = SkyCoord(ra+' '+dec, unit=(u.hourangle, u.deg))
		ras[i] = c.ra.to('deg').value
		decs[i] = c.dec.to('deg').value
	return ras, decs

def crossmatch(sdss, ras = ras, decs=decs):
	inds = []
	for i in range(len(sdss)):	
		inds.append(np.argmin(abs(sdss['ra'][i] - ras) + abs(sdss['dec'][i] - decs)))
	return inds