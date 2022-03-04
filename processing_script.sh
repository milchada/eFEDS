"""REMEMBER TO CHECK THE ENERGY PARAMS, FILE NAMES AND TEL IDs BEFORE RUNNING"""

erosita
evtool eventfiles="fm00_300007_020_EventList_c001.fits" outfile="img1.fits" image=yes emin=0.5 emax=2.0 clobber=yes size='8000'
evtool eventfiles="fm00_300008_020_EventList_c001.fits" outfile="img2.fits" image=yes emin=0.5 emax=2.0 clobber=yes size='8000'
evtool eventfiles="fm00_300009_020_EventList_c001.fits" outfile="img3.fits" image=yes emin=0.5 emax=2.0 clobber=yes size='8000'
evtool eventfiles="fm00_300010_020_EventList_c001.fits" outfile="img4.fits" image=yes emin=0.5 emax=2.0 clobber=yes size='8000'

dmcopy "img1.fits[bin x=-800000:800000:80,y=-800000:800000:80]" img1_cut.fits clobber=yes
dmcopy "img2.fits[bin x=-800000:800000:80,y=-800000:800000:80]" img2_cut.fits clobber=yes
dmcopy "img3.fits[bin x=-800000:800000:80,y=-800000:800000:80]" img3_cut.fits clobber=yes
dmcopy "img4.fits[bin x=-800000:800000:80,y=-800000:800000:80]" img4_cut.fits clobber=yes

#punlearn reproject_image #or delete the parameter file in your home directory (e.g. /home/abogdan/cxcds_param4) if it gives you trouble

reproject_image img1_cut.fits img2_cut.fits img1.fits clobber=yes
reproject_image img3_cut.fits img2_cut.fits img3.fits clobber=yes
reproject_image img4_cut.fits img2_cut.fits img4.fits clobber=yes
mv img2_cut.fits img2.fits

dmimgcalc img1.fits,img2.fits,img3.fits,img4.fits none img_0520.fits op="imgout=img1+img2+img3+img4"

rm img1.fits img2.fits img3.fits img4.fits

expmap inputdatasets="fm00_300007_020_EventList_c001.fits" templateimage="../img_0520_excl_grp_pts.fits" emin=0.5 emax=2.0 mergedmaps='expmap1.fits'
expmap inputdatasets="fm00_300008_020_EventList_c001.fits" templateimage="../img_0520_excl_grp_pts.fits" emin=0.5 emax=2.0 mergedmaps='expmap2.fits'
expmap inputdatasets="fm00_300009_020_EventList_c001.fits" templateimage="../img_0520_excl_grp_pts.fits" emin=0.5 emax=2.0 mergedmaps='expmap3.fits'
expmap inputdatasets="fm00_300010_020_EventList_c001.fits" templateimage="../img_0520_excl_grp_pts.fits" emin=0.5 emax=2.0 mergedmaps='expmap4.fits'

dmimgcalc expmap1.fits,expmap2.fits,expmap3.fits,expmap4.fits none exp_0520.fits op="imgout=img1+img2" clobber=yes

dmcopy "img_0520.fits[exclude pos=region(grp_cat.reg)]" img_0520_12346_excl_grp.fits
dmcopy "img_0520_12346_excl_grp.fits[pos=region(invert_ptsrc_0520.fits)]" img_0520_12346_excl_grp_pts.fits clobber=yes

# dmlist merged_img_excl_prof.fits cols

# dmlist merged_img_excl_prof.fits'[cols counts,cel_area,cel_bri,net_err]' data > put_in_a_file.txt