
erosita
evtool eventfiles="fm00_300007_020_EventList_c001.fits" outfile="img1.fits" image=yes emin=0.5 emax=2.0 clobber=yes size='8000' telid='1 2 3 4 6'
evtool eventfiles="fm00_300008_020_EventList_c001.fits" outfile="img2.fits" image=yes emin=0.5 emax=2.0 clobber=yes size='8000' telid='1 2 3 4 6'
evtool eventfiles="fm00_300009_020_EventList_c001.fits" outfile="img3.fits" image=yes emin=0.5 emax=2.0 clobber=yes size='8000' telid='1 2 3 4 6'
evtool eventfiles="fm00_300010_020_EventList_c001.fits" outfile="img4.fits" image=yes emin=0.5 emax=2.0 clobber=yes size='8000' telid='1 2 3 4 6'

dmcopy "img1.fits[bin x=-800000:800000:80,y=-800000:800000:80]" img1_cut.fits clobber=yes
dmcopy "img2.fits[bin x=-800000:800000:80,y=-800000:800000:80]" img2_cut.fits clobber=yes
dmcopy "img3.fits[bin x=-800000:800000:80,y=-800000:800000:80]" img3_cut.fits clobber=yes
dmcopy "img4.fits[bin x=-800000:800000:80,y=-800000:800000:80]" img4_cut.fits clobber=yes

#punlearn reproject_image #or delete the parameter file in your home directory (e.g. /home/abogdan/cxcds_param4) if it gives you trouble

reproject_image img1_cut.fits exp_0520_12346_excl_grp_pts.fits img1.fits clobber=yes
reproject_image img2_cut.fits exp_0520_12346_excl_grp_pts.fits img2.fits clobber=yes
reproject_image img3_cut.fits exp_0520_12346_excl_grp_pts.fits img3.fits clobber=yes
reproject_image img4_cut.fits exp_0520_12346_excl_grp_pts.fits img4.fits clobber=yes

dmimgcalc img1.fits,img2.fits,img3.fits,img4.fits none img_0520_12346.fits op="imgout=img1+img2+img3+img4"

rm img1.fits img2.fits img3.fits img4.fits

expmap inputdatasets="fm00_300007_020_EventList_c001.fits" templateimage="img_0510_12346.fits" emin=0.5 emax=1.0 withsinglemaps=yes withmergedmaps=no singlemaps='exp11.fits exp12.fits exp13.fits exp14.fits exp15.fits exp16.fits exp17.fits'
expmap inputdatasets="fm00_300008_020_EventList_c001.fits" templateimage="img_0510_12346.fits" emin=0.5 emax=1.0 withsinglemaps=yes withmergedmaps=no singlemaps='exp21.fits exp22.fits exp23.fits exp24.fits exp25.fits exp26.fits exp27.fits'
expmap inputdatasets="fm00_300009_020_EventList_c001.fits" templateimage="img_0510_12346.fits" emin=0.5 emax=1.0 withsinglemaps=yes withmergedmaps=no singlemaps='exp31.fits exp32.fits exp33.fits exp34.fits exp35.fits exp36.fits exp37.fits'
expmap inputdatasets="fm00_300010_020_EventList_c001.fits" templateimage="img_0510_12346.fits" emin=0.5 emax=1.0 withsinglemaps=yes withmergedmaps=no singlemaps='exp41.fits exp42.fits exp43.fits exp44.fits exp45.fits exp46.fits exp47.fits'

#now here I have to be careful. I want to add maps 1, 2, 3 and 4 for T1, T2, T3, T4 and T6; separately for T5 and T7
dmimgcalc exp11.fits,exp12.fits,exp13.fits,exp14.fits,exp16.fits none exp_0510_1.fits op="imgout=img1+img2" clobber=yes
dmimgcalc exp21.fits,exp22.fits,exp23.fits,exp24.fits,exp26.fits none exp_0510_2.fits op="imgout=img1+img2" clobber=yes
dmimgcalc exp31.fits,exp32.fits,exp33.fits,exp34.fits,exp36.fits none exp_0510_3.fits op="imgout=img1+img2" clobber=yes
dmimgcalc exp41.fits,exp42.fits,exp43.fits,exp44.fits,exp46.fits none exp_0510_4.fits op="imgout=img1+img2" clobber=yes
rm exp1*fits exp2*fits exp3*fits exp4*fits

dmimgcalc exp_0510_1.fits,exp_0510_2.fits,exp_0510_3.fits,exp_0510_4.fits none exp_0510_12346.fits op="imgout=img1+img2+img3+img4" clobber=yes
rm exp_0520_1.fits exp_0520_2.fits exp_0520_3.fits exp_0520_4.fits

#confirming that what i ran is:
#expmap inputdatasets="fm00_300010_020_EventList_c001.fits" templateimage="../img_0520_excl_grp_pts.fits" emin=0.5 emax=2.0 withsinglemaps=yes 
								   #withmergedmaps=no singlemaps='exp41.fits exp42.fits exp43.fits exp44.fits exp45.fits exp46.fits exp47.fits'