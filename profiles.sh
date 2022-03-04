#!/bin/bash

set -e
punlearn dmextract
galid=1
while mapfile -t -n 10 gal && ((${#gal[@]})); do
    punlearn dmextract
    input_img="../img_0520_excl_grp_pts.fits[bin sky=@${gal[@]}]"
    output_img=boss_profile_${galid}.fits
    dmextract infile=${input_img} outfile=${output_img} clob+
    echo ${galid} "done"
    let galid=${galid}+1
done < boss_gals.reg