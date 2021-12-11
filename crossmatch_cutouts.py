#this happened because for some reason there are 1344 low-mass cutouts rather than 1233 
#and I can't figure out why
from astropy.coordinates import SkyCoord

lowfiles = glob.glob('lowmass/cutouts/0520/*fits'); lowfiles.sort()
highfiles = glob.glob('highmass/cutouts/0520/*fits'); highfiles.sort()
	
def crossmatch(files, sdss):
	radec_sdss = {}
	for i in range(len(sdss)):
		c = SkyCoord(ra = sdss['ra'][i]*u.degree, dec = sdss['dec'][i]*u.degree)
		c = c.to_string('hmsdms').replace('h',':').replace('d',':').replace('m',':').replace('s ',',').replace('s','')
		radec_sdss[i] = c

	radec = {}
	for f in files:
		num = int(f.split('_')[-1].split('.')[0])
		h = fits.getheader(f)['HISTORY']
		for i in range(len(h)):
			if 'box' in h[i]:
				ra_dec = h[i].split("box(")[1].split('A')[0] + h[i+1].split(' :')[1].split(",1.66'")[0]
				# print(num, ra_dec)
				radec[num] = ra_dec

	ids = {}
	for num in radec_sdss.keys():
		for i in radec.keys():
			k = 1
			ra, dec = radec_sdss[num].split(',')
			if (ra[:6] in radec[i]) and (dec[:6] in radec[i]): #those integers set the precision of the matching
				print(radec[i], radec_sdss[num])
				ids[num] = i
				k = 0
		if k:
			ids[num] = None
	return ids 

lowids = crossmatch(lowfiles, sdss1)
highids = crossmatch(highfiles, sdss2)

#next, I find the cutouts in a mass bin

def get_file(files, num):
	for file in files:
		if '_'+str(num)+'.fits' in file:
			return file

Mstar = np.arange(10.25,11.2,0.1)
dM = 0.05
fig, ax = plt.subplots(nrows=3,ncols=3, sharex=True, sharey=True)
for i in range(5):
	mask  = ((sdss1['logMass'] > (Mstar[i] - dM))*(sdss1['logMass'] < (Mstar[i] + dM)))
	filesub = [get_file(lowfiles, lowids[j]) for j in range(len(mask)) if mask[j]]
	im = stacked_image(filesub, 26, 26)
	img = ax.flatten()[i].imshow(im, cmap=cm.afmhot, origin='lower', norm=colors.Normalize(0,8))
	ax.flatten()[i].set_title('%0.1f' % Mstar[i])
	fig.colorbar(img, ax=ax.flatten()[i])

for i in range(5,len(Mstar)-1):
	mask  = ((sdss2['logMass'] > (Mstar[i] - dM))*(sdss2['logMass'] < (Mstar[i] + dM)))
	filesub = [get_file(highfiles,highids[j]) for j in range(len(mask)) if mask[j]]
	im = stacked_image(filesub, 26, 26)
	img = ax.flatten()[i].imshow(im, cmap=cm.afmhot, origin='lower', norm=colors.Normalize(0,8))
	ax.flatten()[i].set_title('%0.1f' % Mstar[i])
	fig.colorbar(img, ax=ax.flatten()[i])