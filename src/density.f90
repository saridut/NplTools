program density
    !! Calculates the axial distribution function of specified atom types from
    !! a LAMMPS binary trajectory file.

use M_CLI2, only : set_args, get_args, get_args_fixed_length, &
                    get_args_fixed_size
use constants_m
use strings_m
use trajectory_reader_m
use simbox_m
use atmcfg_m
use config_io_m
use zdf_calculator_m

implicit none

character(len=80), allocatable, dimension(:) :: help_text
character(2) :: coord_desc(3), boundary(3)
character(:), allocatable :: fn_top, fn_molinfo, fn_fns, mol_names(:)

real(rp) :: box(9)
real(rp), allocatable :: dmp_data(:,:)

integer(ip_long) :: nts, natoms
integer :: col_nums(5), ifrm_beg, ifrm_end, ifrm_stp, i, iatm, ierr, &
            iframe, atm_id, at
integer, allocatable :: mol_bnds(:,:)

logical :: is_triclinic, read_markers

type(smbx_t) :: simbox
type(atmcfg_t) :: atc
type(trajectory_reader_t) :: trdr
type(zdf_calculator_t)    :: zdfc

call set_args('         &
    & --fn_top " "      &
    & --fn_fns " "      &
    & --read_markers  F &
    & --fn_molinfo " "  &
    & --col_nums []     &
    & --coord_desc "sw","sw","sw"  &
    & --fn_out " "      &
    & --frame 1         &
    & --to -1           &
    & --stp 1           &
    & --is_slab F       &
    & --normalize F     &
    & --mass_density F  &
    & --zmax 0.0        &
    & --zads 0.0        &
    & --bin_size 0.0    &
    & --ats_xtal []     &
    & ', help_text)


call get_args('fn_top', fn_top)
call get_args('fn_fns', fn_fns)
call get_args('read_markers', read_markers)
call get_args('fn_molinfo', fn_molinfo)

call get_args_fixed_size('col_nums', col_nums)
call get_args_fixed_size('coord_desc', coord_desc)

call get_args('fn_out', fn_out)

call get_args('frame', ifrm_beg)
call get_args('to', ifrm_end)
call get_args('stp', ifrm_stp)

write(*,"(a,1x,a)") 'fn_top =', fn_top
write(*,"(a,1x,a)") 'fn_fns =', fn_fns
write(*,"(a,1x,l1)") 'read_markers =', read_markers
write(*,"(a,1x,a)") 'fn_molinfo =', fn_molinfo
write(*,"(a,1x,*(i0,1x))") 'col_nums =', col_nums
write(*,"(a,1x,*(a,1x))") 'coord_desc =', coord_desc
write(*,"(a,1x,a)") 'fn_out =', fn_out
write(*,"(3(a,i0,1x),/)") 'frame = ', ifrm_beg, 'to = ', ifrm_end, &
    'stp = ', ifrm_stp

call simbox%init()
call read_ldf(simbox, atc, fn_top)
atc%coord_desc = coord_desc
call read_molinfo(fn_molinfo, mol_names, mol_bnds)

call zdfc%init(atc)

!Open trajectory reader
call trdr%init(fn_fns, col_nums, ierr, read_markers)
if (ierr /= 0) error stop
if (.not. read_markers) call trdr%write_markers(fn_fns)

write(*,'(a,i0)') 'numframes = ', trdr%num_frames
if (ifrm_beg == -1) ifrm_beg = trdr%num_frames
if (ifrm_end == -1) ifrm_end = trdr%num_frames

do iframe = ifrm_beg, ifrm_end, ifrm_stp
    !Load a frame
    call trdr%read_frame(iframe, nts, natoms, is_triclinic, boundary, box, &
        dmp_data, ierr)
    if (ierr /= 0) exit

    call simbox%set_boundary(boundary(1), boundary(2), boundary(3))
    call simbox%set_bounds(xlo=box(1), xhi=box(2), ylo=box(3), yhi=box(4), &
        zlo=box(5), zhi=box(6))
    if (is_triclinic) then
        simbox%is_triclinic = .true.
        call simbox%set_tilt(xy=box(7), xz=box(8), yz=box(9))
    end if
    call simbox%update()

    do iatm = 1, int(natoms,ip)
        atm_id = int(dmp_data(1,iatm)); at = int(dmp_data(2,iatm))
        if (atm_id <= atc%num_atoms) then
            atc%atoms(1,atm_id) = at
            atc%coordinates(:,atm_id) = dmp_data(3:5,iatm)
        end if
    end do
    call simbox%unscale_all(atc%coordinates, atc%coord_desc)

    call zdfc%process_frame(atc, mol_names, mol_bnds)

    if (mod((iframe-ifrm_beg),10)==0) write(*,"(a,1x,i0,1x,a)") "Frame", iframe, "done"
end do

!Close trajectory reader
call trdr%delete()

!Write out data
call zdfc%writeout()
call zdfc%delete()

!*******************************************************************************

end program
