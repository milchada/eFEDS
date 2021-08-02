#!/bin/bash

set -e

galid=1
cat highmass_boxes.reg | while read -r gal; 
   do
	  input_img='../merged_img_excl_grp_pts.fits'
	  output_img=highmass_cutout_${galid}.fits
          
	  dmcopy infile="${input_img}[pos=${gal}]" \
		  outfile=${output_img} clob+
          echo ${galid} "done"
	  let galid=${galid}+1
   done
   