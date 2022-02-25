def Lx_stack(Mstar, cts, sdss, area, condn=None, exp=None, bmin=1,bckcols=3, 
	ret=True, simsfile=None, suffix='', xrb=False):
	dM = (Mstar[1] - Mstar[0])/2.
	ms = sdss['logMass']
	Lx = np.zeros(len(Mstar))
	binerr = np.zeros(len(Mstar))

	for i in range(len(Mstar)):
		mask  = (ms > (Mstar[i] - dM))*(ms < (Mstar[i] + dM))
		mask = [j for j in range(len(mask)) if mask[j]]

		if sum(mask):
			sfr = sdss[mask]['ssfr'] + ms[mask]
			if condn:
				if (condn == 'up2'):
					cond = (sfr > np.nanmedian(sfr))
				elif (condn=='down2'):
					cond = (sfr < np.nanmedian(sfr))
				elif (condn == 'up3') or (condn=='down3'):
					sort = np.argsort(sfr)
					ind1 = int(len(sfr)/3)
					print(sfr[ind1])
					if condn == 'up3':
						cond = (sfr > sfr[sort[int(2*ind1)]])
					else:
						cond = (sfr < sfr[sort[ind1]])
				scts = cts[mask][cond]
				sexp = exp[mask][cond]
				sz = sdss[mask][cond]
				sarea = area[mask][cond]
			else:
				scts = cts[mask]
				sexp = exp[mask]
				sz = sdss[mask]
				sarea = area[mask]
			print(len(scts))
			if simsfile:
				fig, ax = profiles(scts, sexp, sarea, sz,bckcols=bckcols, ret=False, simsfile=simsfile, xrb=xrb)
				fig.savefig('profile_%0.1f_%0.1f%s.png' % (Mstar[i]-dM, Mstar[i]+dM, suffix))
			else:
				x, sb, binerr[i], sb_bkg = profiles(scts, sexp, sarea,sz, bckcols=bckcols, simsfile=simsfile, xrb=xrb)
				print(sb, sb_bkg)
				x1 = x+10
				ar = np.pi*(x1**2 - x**2)
				ypos = (sb*ar)#[sb > 0.05*sb_bkg] #select the bins where signal > 0.05*err
				print(ypos.shape)

				Lx[i] = sum(ypos[bmin:-bckcols]) #exclude central bin

	if ret:
		return Lx, binerr

def Lx_obs(cts, area, exp, sdss, split = False, dM = 0.1, bckcols=5, xrb=True):	
	Mstar = np.arange(10.2+dM,11.2,2*dM)
	
	if split in [3, 4]:
		sfr = 10**(sdss['ssfr'] + sdss['logMass'])
	if split == 1:
		Lx_blue, binerr_blue = Lx_stack(Mstar, cts, sdss, area, condn = 'up2', exp=exp, xrb=xrb) 
		Lx_red, binerr_red = Lx_stack(Mstar, cts, sdss, area, condn = 'down2', exp=exp, xrb=xrb) 
	elif split == 2:	
		Lx_blue, binerr_blue = Lx_stack(Mstar, cts, sdss, area, condn = 'up3', exp=exp, xrb=xrb) 
		Lx_red, binerr_red = Lx_stack(Mstar, cts, sdss, area, condn = 'down3', exp=exp, xrb=xrb) 
	elif split == 3:	
		Lx_blue, binerr_blue = Lx_stack(Mstar, cts[sfr>0.1], sdss[sfr>0.1], area[sfr>0.1], exp=exp[sfr>0.1], xrb=xrb) 
		Lx_red, binerr_red = Lx_stack(Mstar, cts[sfr<0.1], sdss[sfr<0.1], area[sfr<0.1], exp=exp[sfr<0.1], xrb=xrb) 
	elif split == 4:
		Lx_blue, binerr_blue = Lx_stack(Mstar, cts[sfr>1], sdss[sfr>1], area[sfr>1], exp=exp[sfr>1], xrb=xrb) 
		Lx_red, binerr_red = Lx_stack(Mstar, cts[sfr<1], sdss[sfr<1.0], area[sfr<1.0], exp=exp[sfr<1.0], xrb=xrb) 
	if split:
		return Lx_blue, binerr_blue, Lx_red, binerr_red 
	else:
		Lx, binerr = Lx_stack(Mstar, cts, sdss, area, condn=None, exp=exp, xrb=xrb) 
		return Lx, binerr
	