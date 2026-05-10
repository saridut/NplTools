#!/usr/bin/env bash

#fn_cfg="/Volumes/lasdata/ML3_20x20_ac/NPL_ML3_20x20_45_ac.lmp"
#fn_traj="/Volumes/lasdata/ML3_20x20_ac/trajeq_lig_0.bin"

#fn_cfg="/Volumes/lasdata/ML3_20x20_ol/NPL_ML3_20x20_45_ol.lmp"
#fn_traj="/Volumes/lasdata/ML3_20x20_ol/trajeq_lig_0.bin"

fn_cfg="/Volumes/lasdata/ML3_20x20/NPL_ML3_20x20_45.lmp"
fn_traj="/Volumes/lasdata/ML3_20x20/trajeq_lig_0.bin"

#fn_mol="molrec_ac.txt"
#fn_mol="molrec_ol.txt"
fn_mol="molrec_olac.txt"

#Number of the first frame
ifrm_beg=1
#Number of the last frame. -1 indicates the last available frame.
ifrm_end=-1
#Step over how many frames? 1 indicates consecutive frames.
ifrm_stp=1

r_ads='2.0'
bs_theta='2'
bs_phi='2'

./calc_angle $fn_cfg $fn_mol $fn_traj $ifrm_beg $ifrm_end $ifrm_stp $r_ads \
    $bs_theta $bs_phi
