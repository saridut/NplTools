program calc_angle
    !! Calculates the binding motif (angles) for adsorbing ligands from
    !! a LAMMPS binary trajectory file.

use constants_m
use strings_m
use utils_math_m
use trajectory_m
use simbox_m
use atmcfg_m
use config_io_m
!use cell_list_m

implicit none

character(len=256) :: cla_buf
character(len=:), allocatable :: fn_cfg, fn_mol, fn_traj

type(smbx_t) :: simbox
type(atmcfg_t) :: atc
type(trajectory_t) :: traj

real(rp), dimension(:), allocatable :: buf_pad
real(rp), dimension(9) :: box
integer, dimension(2,3) :: boundary
integer :: buf_pad_len, atm_id, at, num_fields
integer(ip_long) :: nts

integer :: ifrm_beg, ifrm_end, ifrm_stp
integer, dimension(:,:), allocatable :: molecules

real(rp), dimension(:), allocatable :: bins_phi, bins_theta
real(rp), dimension(:), allocatable :: hist_phi, hist_theta
real(rp), dimension(2,3) :: xtal_xtnt ! Extent of the crystal
real(rp), dimension(3) :: ri, rj, rij, r_C, r_O1, r_O2, r_C_O1, r_C_O2, uhat
real(rp) :: rij_mag, zdist !For radial and axial distances
real(rp) :: r_ads ! Distance cutoff for adsorption
real(rp) :: theta, phi, binsiz_theta, binsiz_phi
integer :: num_molecules, num_atoms_xtal, num_bins_theta, num_bins_phi
integer :: iatm, imol, ierr, ibeg, mt, iatmC
integer :: iframe, oframe, nofrms
integer :: fu_mol, fu_out
integer :: i, j, k

call get_command_argument(1, cla_buf)
fn_cfg = trim(cla_buf)

call get_command_argument(2, cla_buf)
fn_mol = trim(cla_buf)

call get_command_argument(3, cla_buf)
fn_traj = trim(cla_buf)

call get_command_argument(4, cla_buf)
ifrm_beg = int(str_to_d(cla_buf))

call get_command_argument(5, cla_buf)
ifrm_end = int(str_to_d(cla_buf))

call get_command_argument(6, cla_buf)
ifrm_stp = int(str_to_d(cla_buf))

call get_command_argument(7, cla_buf)
r_ads = str_to_d(cla_buf)

call get_command_argument(8, cla_buf)
binsiz_theta = (math_pi/180.0_rp)*str_to_d(cla_buf)

call get_command_argument(9, cla_buf)
binsiz_phi = (math_pi/180.0_rp)*str_to_d(cla_buf)

!Read in LAMMPS data file containing the initial configuration
call read_ldf(simbox, atc, fn_cfg)

!Read the molecular data
open(newunit=fu_mol, file=fn_mol, action='read', status='old')
read(fu_mol,*) num_molecules
allocate(molecules(3,num_molecules))
read(fu_mol,*)
do i = 1, num_molecules
    read(fu_mol,*) imol, molecules(:,i)
end do
close(fu_mol)

!Determine the extent of the crystal
xtal_xtnt(1,:) = huge(1.0_rp); xtal_xtnt(2,:) = -huge(1.0_rp)
num_atoms_xtal = 0
do iatm = 1, atc%num_atoms
    if (atc%atoms(1,iatm) < 3) then
        num_atoms_xtal = num_atoms_xtal + 1
        ri = atc%coordinates(:,iatm)
        do i = 1, 3
            xtal_xtnt(1,i) = min( xtal_xtnt(1,i), ri(i) )
            xtal_xtnt(2,i) = max( xtal_xtnt(2,i), ri(i) )
        end do
    end if
end do

!Set up bins for histograms: theta
num_bins_theta = ceiling(math_pi_2/binsiz_theta)
binsiz_theta = math_pi_2/num_bins_theta

allocate(bins_theta(num_bins_theta))
allocate(hist_theta(num_bins_theta))
hist_theta = 0.0_rp

do i = 1, num_bins_theta
    bins_theta(i) = i*binsiz_theta- 0.5*binsiz_theta
end do
!Set up bins for histograms: phi
num_bins_phi = ceiling(math_pi/binsiz_phi)
binsiz_phi = math_pi/num_bins_phi

allocate(bins_phi(num_bins_phi))
allocate(hist_phi(num_bins_phi))
hist_phi = 0.0_rp

do i = 1, num_bins_phi
    bins_phi(i) = -math_pi_2 + i*binsiz_phi - 0.5*binsiz_phi
end do

write(*,'(a,2x,i0)')   'NUM_ATMS_XTAL', num_atoms_xtal
write(*,'(a,2x,g0.8,1x,g0.8)') 'XTAL_X', xtal_xtnt(:,1)
write(*,'(a,2x,g0.8,1x,g0.8)') 'XTAL_Y', xtal_xtnt(:,2)
write(*,'(a,2x,g0.8,1x,g0.8)') 'XTAL_Z', xtal_xtnt(:,3)

!Chop off the edges of the crytal
xtal_xtnt(1,1:2) = xtal_xtnt(1,1:2) + 6.0_rp
xtal_xtnt(2,1:2) = xtal_xtnt(2,1:2) - 6.0_rp

!------------------------------------------------------------------------------
!Open file for reading trajectory
call traj%open(fn_traj, ierr)
write(*,*) 'numframes ', traj%num_frames

buf_pad_len = int( traj%buf_pad_size/sizeof_real )
allocate(buf_pad(buf_pad_len))

call simbox%init(is_bounded=.true., is_triclinic=traj%is_triclinic)

if (ifrm_beg == -1) ifrm_beg = traj%num_frames
if (ifrm_end == -1) ifrm_end = traj%num_frames
nofrms = (ifrm_end - (ifrm_beg-1))/ifrm_stp

oframe = 0
do iframe = ifrm_beg, ifrm_end, ifrm_stp
    call traj%read(iframe, nts, box, boundary, buf_pad, ierr)
    oframe = oframe + 1
    if (mod((iframe-ifrm_beg),10)==0) write(*,*) oframe, iframe

    !Set the box bounds and tilt
    call simbox%set_bound(xlo=box(1), xhi=box(2), &
        ylo=box(3), yhi=box(4), zlo=box(5), zhi=box(6))
    if (traj%is_triclinic) then
        call simbox%set_tilt(xy=box(7), xz=box(8), yz=box(9))
    end if
    call simbox%update()

    !Read in the atom coordinates
    num_fields = traj%size_one
    do iatm = 1, traj%num_atoms
        ibeg = (iatm-1)*num_fields + 1
        atm_id = int(buf_pad(ibeg)); at = int(buf_pad(ibeg+1))
        atc%atoms(1,atm_id) = at
        atc%coordinates(:,atm_id) = buf_pad(ibeg+2:ibeg+4)
    end do

    !Unscale atom coordinates
    call simbox%unscale_all(atc%coordinates(:,num_atoms_xtal+1:atc%num_atoms))

    !Loop over all ligand molecules
    do imol = 1, num_molecules
        mt = molecules(1,imol); ibeg = molecules(2,imol)
        if (mt == 1) iatmC = ibeg
        if (mt == 2) iatmC = ibeg + 17

        !Check if the molecule is adsorbed
        r_C = atc%coordinates(:,iatmC)

        if ( (r_C(1) < xtal_xtnt(1,1)) .or. (r_C(1) > xtal_xtnt(2,1)) ) cycle
        if ( (r_C(2) < xtal_xtnt(1,2)) .or. (r_C(2) > xtal_xtnt(2,2)) ) cycle

        if ( r_C(3) >= xtal_xtnt(2,3) ) then
            zdist = r_C(3) - xtal_xtnt(2,3)
        else if ( r_C(3) <= xtal_xtnt(1,3) ) then
            zdist = xtal_xtnt(1,3) - r_C(3)
        end if
        if (zdist > r_ads) cycle

        !Get the positions of the two oxygen atoms
        if (mt == 1) then
            r_O1 = atc%coordinates(:,ibeg+2)
            r_O2 = atc%coordinates(:,ibeg+3)
        else if (mt == 2) then
            r_O1 = atc%coordinates(:,ibeg+18)
            r_O2 = atc%coordinates(:,ibeg+19)
        end if
        r_C_O1 = r_O1 - r_C
        r_C_O2 = r_O2 - r_C
        call cross(r_C_O1, r_C_O2, uhat)
        uhat = uhat/norm2(uhat)
        theta = acos(uhat(3)); phi = atan(uhat(2)/uhat(1))
        if (theta > math_pi_2) theta = math_pi - theta
        k = floor(theta/binsiz_theta) + 1; hist_theta(k) = hist_theta(k) + 1
        k = floor((phi+math_pi_2)/binsiz_phi) + 1; hist_phi(k) = hist_phi(k) + 1
    end do !End of loop over ligands

end do !End of loop over frames

hist_theta = hist_theta/nofrms !Average over all frames
if (any(hist_theta>0)) hist_theta = hist_theta/(binsiz_theta*sum(hist_theta)) !Normalize

hist_phi = hist_phi/nofrms !Average over all frames
if (any(hist_phi>0)) hist_phi = hist_phi/(binsiz_phi*sum(hist_phi)) !Normalize

!Close trajectory
call traj%close()

!Writing out histogram data
open(newunit=fu_out, file='hist_theta.csv', action='write', status='replace')
write(fu_out, '(a16,",",1x,a16)') 'theta', 'hist'

bins_theta = bins_theta*180._rp/math_pi
do i = 1, num_bins_theta
    write(fu_out, '(es16.7,",",1x,es16.7)') bins_theta(i), hist_theta(i)
end do
close(fu_out)

bins_phi = bins_phi*180._rp/math_pi
open(newunit=fu_out, file='hist_phi.csv', action='write', status='replace')
write(fu_out, '(a16,",",1x,a16)') 'phi', 'hist'

do i = 1, num_bins_phi
    write(fu_out, '(es16.7,",",1x,es16.7)') bins_phi(i), hist_phi(i)
end do
close(fu_out)

!*******************************************************************************

end program
