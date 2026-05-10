program bin2txt

use M_CLI2, only : set_args, get_args, get_args_fixed_length, &
                    get_args_fixed_size
use constants_m
use strings_m
use trajectory_reader_m
use simbox_m
use atmcfg_m
use config_io_m

implicit none

character(len=80), allocatable, dimension(:) :: help_text
character(len=:), allocatable :: fn_cfg
    !! Name of LAMMPS data file containing a reference configuration with topology
character(len=:), allocatable :: fn_fns
    !! Name of file containing the list of dump file names
character(len=:), allocatable :: fn_markers
    !! Name of file to write out frame markers
character(len=:), allocatable :: fn_out
    !! Output filename pattern containing one and only one `*` which will be 
    !! replaced by the value of the corresponding time step.
character(len=2), dimension(3) :: coord_desc
character(len=2), dimension(:), allocatable :: atm_names
integer :: ifrm_beg, ifrm_end, ifrm_stp

type(smbx_t) :: simbox
type(atmcfg_t) :: atc
type(trajectory_reader_t) :: trdr
logical :: is_triclinic, read_markers, write_markers
integer(ip_long) :: nts, natoms
character(len=2), dimension(3) :: boundary
integer, dimension(5) :: col_nums
real(rp), dimension(9) :: box
real(rp), dimension(:,:), allocatable :: dmp_data
character(len=3) :: file_type, clinp
character(len=:), allocatable :: fn
integer :: i, iatm, ierr, ipos, iframe, oframe, nofrms, atm_id, at

help_text = [&
'USAGE                                                                       ',&
'  bin2txt --fn_top <file> --fn_fns <file> --col_nums <list of ints>         ',&
'    --coord_desc <list of strings> --atm_names <list of strings>            ',&
'    --fn_out <string> --frame <int> --to <int> --stp <int>                  ',&
'    [--read_markers] [--write_markers] [--fn_markers <file>]                ',&
'                                                                            ',&
'DESCRIPTION                                                                 ',&
'  Converts LAMMPS binary dump file to a single LAMMPS ascii dump file, or a ',&
'  series of files (each per time step) in LAMMPS ascii dump format, LAMMPS  ',&
'  data file format, or XYZ format.                                          ',&
'                                                                            ',&
'OPTIONS                                                                     ',&
'  --fn_top <file>                                                           ',&
'    File containing the intial topology in LAMMPS data file format.         ',&
'                                                                            ',&
'  --fn_fns <file>                                                           ',&
'    File containing the names of the trajectory files.                      ',&
'                                                                            ',&
'  --col_nums <list of ints>                                                 ',&
'    Comma-separated list of five or eight integers specifying the index     ',&
'    (1-based) of the following quantities in the per-atom record of the dump',&
'    file: atom id, atom type, position(x), position(y), position(z),        ',&
'    image(x), image(y), image(z). The last three may be omitted.            ',&
'                                                                            ',&
'  --coord_desc <list of strings> [default: "sw","sw","sw"]                  ',&
'    Comma-separated list of three 2-character strings representing scaling  ',&
'    and wrapping of the atom coordinates along the x, y, and z directions,  ',&
'    respectively. The first character of each string must be either `u`     ',&
'    (unscaled) or `s` (scaled) and the second character must be either `u`  ',&
'    (unwrapped) or `w`(wrapped).                                            ',&
'                                                                            ',&
'  --atm_names <list of strings>                                             ',&
'    Names of all atoms in the topolgy file as a comma-separated list of     ',&
'    2-character strings. Must be present if converting to XYZ format, may be',&
'    omitted for other formats.                                              ',&
'                                                                            ',&
'  --fn_out <string>                                                         ',&
'     Output file name with extension `.txt` for LAMMPS ascii dump format,   ',&
'     `.lmp` for LAMMPS data file format, and `.xyz` for XYZ format. To have ',& 
'     a separate file for each time step, the file name must contain at most ',&
'     one `*`, which will be replaced by the corresponding time step number. ',&
'                                                                            ',&
'  --frame <int> [default: 1]                                                ',&
'    Frame number to begin with (1-based).                                   ',&
'                                                                            ',&
'  --to <int> [default: -1]                                                  ',&
'    Frame number to end at (1-based). -1 indicates the last frame of the    ',&
'    trajectory.                                                             ',&
'                                                                            ',&
'  --stp <int> [default: 1]                                                  ',&
'    Frame step size.                                                        ' &
    ]

call set_args('         &
    & --fn_top " "      &
    & --fn_fns " "      &
    & --col_nums []     &
    & --coord_desc "sw","sw","sw"  &
    & --atm_names []    &
    & --fn_out " "      &
    & --frame 1         &
    & --to -1           &
    & --stp 1           &
    & --read_markers  F &
    & --write_markers F &
    & --fn_markers " "  &
    & ', help_text)


call get_args('fn_top', fn_top)
call get_args('fn_fns', fn_fns)

call get_args_fixed_size('col_nums', col_nums)
call get_args_fixed_size('coord_desc', coord_desc)

call get_args_fixed_length('atm_names', atm_names)
call get_args('fn_out', fn_out)

call get_args('frame', ifrm_beg)
call get_args('to', ifrm_end)
call get_args('stp', ifrm_stp)

call get_args('read_markers', read_markers)
call get_args('write_markers', write_markers)
call get_args('fn_markers', fn_markers)

write(*,"(a,1x,a)") 'fn_top =', fn_top
write(*,"(a,1x,a)") 'fn_fns =', fn_fns
write(*,"(a,1x,*(i0,1x))") 'col_nums[5] =', col_nums
write(*,"(a,1x,*(a,1x))") 'coord_desc[3] =', coord_desc

write(*,"(a,i0,a,*(a,1x))") 'atm_names[', size(atm_names), '] = ', &
    (trim(atm_names(i)), i=1, size(atm_names))

write(*,"(a,1x,a)") 'fn_out =', fn_out
write(*,"(3(a,i0,1x))") 'ifrm_beg = ', ifrm_beg, 'ifrm_end = ', ifrm_end, &
    'ifrm_stp = ', ifrm_stp
write(*,"(a,1x,l1)") 'read_markers =', read_markers
write(*,"(a,1x,l1)") 'write_markers =', write_markers
write(*,"(a,1x,a)") 'fn_markers =', fn_markers

if (str_endswith(fn_out, '.xyz')) then
    file_type = 'xyz'
else if (str_endswith(fn_out, '.lmp')) then
    file_type = 'ldf'
else if (str_endswith(fn_out, '.txt')) then
    file_type = 'dmp'
else
    write(*,"(a,1x,a)") "ERROR: Unrecognized file extension for", fn_out
    error stop
end if

ipos = scan(fn_out, '*')


call simbox%init()
call read_ldf(simbox, atc, fn_top)
atc%coord_desc = coord_desc

if ( allocated(atm_names) .and. (size(atm_names)>0) ) then
    do i = 1, size(atm_names)
        atc%atom_specs(i)%name = trim(atm_names(i))
    end do
end if

!Open trajectory reader.
call trdr%init(fn_fns, col_nums, ierr, read_markers)
if (ierr /= 0) error stop

if (write_markers) call trdr%write_markers(fn_markers)

write(*,'(a,i0)') 'numframes = ', trdr%num_frames

if (ifrm_beg == -1) ifrm_beg = trdr%num_frames
if (ifrm_end == -1) ifrm_end = trdr%num_frames

!Number of output frames
nofrms = (ifrm_end - (ifrm_beg-1))/ifrm_stp
if (nofrms > 1000) then
    write(*,*) 'More than 1000 frames will be written'
    write(*,'(a)', advance='no') 'Press "y" if ok: '
    read(*,*) clinp
    !Will fall through if not 'y' 
    if ( trim(clinp) /= 'y' ) ifrm_end = ifrm_beg - 1
end if

oframe = 0
do iframe = ifrm_beg, ifrm_end, ifrm_stp
    call trdr%read_frame(iframe, nts, natoms, is_triclinic, boundary, box, &
        dmp_data, ierr)
    if (ierr /= 0) exit
    oframe = oframe + 1
    if (mod((iframe-ifrm_beg),1)==0) write(*,*) oframe, iframe, nts

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

    call simbox%unscale_all(atc%coordinates, coord_desc)

    fn = fn_out(1:ipos-1)//str_from_num(nts)//fn_out(ipos+1:)
    if (file_type == 'xyz') then
        call write_xyz(atc, fn, 'timestep = '//str_from_num(nts))
    else if (file_type == 'ldf') then
        call write_ldf(simbox, atc, fn, 'timestep = '//str_from_num(nts))
    else if (file_type == 'dmp') then
        call write_dmp(simbox, atc, nts, fn)
    end if
end do

!Close trajectory reader
call trdr%delete()

!*******************************************************************************

end program
