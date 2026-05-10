program calc_rdf
    !! Calculates the radial/axial distribution function of specified atom types from
    !! a LAMMPS binary trajectory file.

use constants_m
use strings_m
use trajectory_m
use simbox_m
use atmcfg_m
use config_io_m

implicit none

character(len=256) :: cla_buf
character(len=:), allocatable :: fn_cfg
character(len=:), allocatable :: fn_traj
character(len=:), allocatable :: fn_out

type(smbx_t) :: simbox
type(atmcfg_t) :: atc
type(trajectory_t) :: traj

real(rp), dimension(:), allocatable :: buf_pad
real(rp), dimension(9) :: box
integer, dimension(2,3) :: boundary

integer :: ifrm_beg, ifrm_end, ifrm_stp
integer(ip_long) :: nts
integer :: buf_pad_len, atm_id, at, num_fields
integer :: iatm, ierr, ibeg
integer :: iframe, oframe, nofrms

real(rp), dimension(:), allocatable :: bins !Bins: (r-dr/2, r+dr/2)
real(rp), dimension(:), allocatable :: gr !Binned RDF
real(rp), dimension(2,3) :: xtal_xtnt ! Extent of the crystal
real(rp), dimension(3) :: ri, rj, rij
real(rp) :: rij_mag, zi_mag !For radial and axial distances
real(rp) :: r_max, dr
integer :: num_bins, num_atoms_xtal
integer :: i, j, k, fu_out

call get_command_argument(1, cla_buf)
fn_cfg = trim(cla_buf)

call get_command_argument(2, cla_buf)
fn_traj = trim(cla_buf)

call get_command_argument(3, cla_buf)
ifrm_beg = int(str_to_d(cla_buf))

call get_command_argument(4, cla_buf)
ifrm_end = int(str_to_d(cla_buf))

call get_command_argument(5, cla_buf)
ifrm_stp = int(str_to_d(cla_buf))

call get_command_argument(6, cla_buf)
r_max = str_to_d(cla_buf)

call get_command_argument(7, cla_buf)
dr = str_to_d(cla_buf)

call get_command_argument(8, cla_buf)
fn_out = trim(cla_buf)

!Read in LAMMPS data file containing the initial configuration
call read_ldf(simbox, atc, fn_cfg)

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

num_bins = ceiling(r_max/dr)
dr = r_max/num_bins !Fix dr to make all bins equispaced

allocate(bins(1:num_bins))
allocate(gr(1:num_bins))

bins = 0.0_rp; gr = 0.0_rp

do i = 1, num_bins
    bins(i) = i*dr - 0.5*dr
end do

write(*,'(a,2x,g0.8)') 'BIN_SIZE', dr
write(*,'(a,2x,g0.8)') 'R_MAX', r_max
write(*,'(a,2x,g0.8)') 'NUM_BINS', num_bins
write(*,'(a,2x,i0)')   'NUM_ATMS_XTAL', num_atoms_xtal
write(*,'(a,2x,g0.8,1x,g0.8)') 'XTAL_X', xtal_xtnt(:,1)
write(*,'(a,2x,g0.8,1x,g0.8)') 'XTAL_Y', xtal_xtnt(:,2)
write(*,'(a,2x,g0.8,1x,g0.8)') 'XTAL_Z', xtal_xtnt(:,3)

!Open file for writing rdf data & write header line
open(newunit=fu_out, file=fn_out, action='write', status='replace')
write(fu_out, '(a16,",",1x,a16)') 'r', 'gr'

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

xtal_xtnt(1,1:2) = xtal_xtnt(1,1:2) + 5.0_rp
xtal_xtnt(2,1:2) = xtal_xtnt(2,1:2) - 5.0_rp

oframe = 0
do iframe = ifrm_beg, ifrm_end, ifrm_stp
    call traj%read(iframe, nts, box, boundary, buf_pad, ierr)
    oframe = oframe + 1
!   if (mod((iframe-ifrm_beg),10)==0) write(*,*) oframe, iframe

    call simbox%set_bound(xlo=box(1), xhi=box(2), &
        ylo=box(3), yhi=box(4), zlo=box(5), zhi=box(6))
    if (traj%is_triclinic) then
        call simbox%set_tilt(xy=box(7), xz=box(8), yz=box(9))
    end if
    call simbox%update()

    num_fields = traj%size_one
    do iatm = 1, traj%num_atoms
        ibeg = (iatm-1)*num_fields + 1
        atm_id = int(buf_pad(ibeg)); at = int(buf_pad(ibeg+1))
        atc%atoms(1,atm_id) = at
        atc%coordinates(:,atm_id) = buf_pad(ibeg+2:ibeg+4)
    end do

    call simbox%unscale_all(atc%coordinates(:,num_atoms_xtal+1:atc%num_atoms))

    !Calculating RDF
!   do iatm = 1, atc%num_atoms-1
!       ri = atc%coordinates(:,iatm)
!       do jatm = iatm+1, atc%num_atoms
!           rj = atc%coordinates(:,jatm)
!           rij = rj - ri
!           rij_mag = sqrt(rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3))
!           if (rij_mag > r_max) cycle
!           k = floor(rij_mag/dr) + 1
!           gr(k) = gr(k) + 2
!       end do
!   end do

    do iatm = num_atoms_xtal+1, atc%num_atoms
        !Only atom type = 3
        at = atc%atoms(1,iatm)
        !if ( (at/=6) .and. (at/=10) ) cycle
        !if ( at/=10 ) cycle
        ri = atc%coordinates(:,iatm)
        if ( (ri(1) < xtal_xtnt(1,1)) .or. (ri(1) > xtal_xtnt(2,1)) ) cycle
        if ( (ri(2) < xtal_xtnt(1,2)) .or. (ri(2) > xtal_xtnt(2,2)) ) cycle

        if ( ri(3) >= xtal_xtnt(2,3) ) then
            zi_mag = ri(3) - xtal_xtnt(2,3)
            !cycle
        else if ( ri(3) <= xtal_xtnt(1,3) ) then
            zi_mag = xtal_xtnt(1,3) - ri(3)
            !cycle
        end if

        if (zi_mag > r_max) cycle
        k = floor(zi_mag/dr) + 1
        gr(k) = gr(k) + 1
    end do

end do

!Close trajectory
call traj%close()

gr = gr/nofrms !Average over all frames
gr = gr/(dr*sum(gr)) !Normalize

do i = 1, num_bins
    !write(fu_out, '(es16.7,",",1x,es16.7)') bins(i), gr(i)
    write(fu_out, '(g0.8,",",1x,g0.8)') bins(i), gr(i)
end do

close(fu_out)

!*******************************************************************************

end program
