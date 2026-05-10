!********************************************************************************!
!                                                                                !
! The MIT License (MIT)                                                          !
!                                                                                !
! Copyright (c) 2024 Sarit Dutta <saridut@gmail.com>                             !
!                                                                                !
! Permission is hereby granted, free of charge, to any person obtaining a copy   !
! of this software and associated documentation files (the "Software"), to deal  !
! in the Software without restriction, including without limitation the rights   !
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      !
! copies of the Software, and to permit persons to whom the Software is          !
! furnished to do so, subject to the following conditions:                       !
!                                                                                !
! The above copyright notice and this permission notice shall be included in all !
! copies or substantial portions of the Software.                                !
!                                                                                !
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     !
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       !
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    !
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         !
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  !
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  !
! SOFTWARE.                                                                      !
!                                                                                !
!********************************************************************************!

module trajectory_reader_m
    !! Routines for reading frames from a trajectory files.
 
use constants_m
use strings_m

implicit none

private
public :: trajectory_reader_t

type trajectory_reader_t
    integer :: num_frames = 0
    integer, dimension(:), allocatable :: col_nums
        !! Column numbers to be read from dump files
    character(len=:), dimension(:), allocatable :: fns_dump
        !! Dump file names
    integer(ip_long), dimension(:,:), allocatable :: frame_markers
        !! (2,*num_frames*) Position corresponding to the start of each frame
    contains
        procedure :: init       => trdr_init
        procedure :: delete     => trdr_delete
        procedure :: read_frame => trdr_read_frame
        procedure :: write_markers => trdr_write_markers
        procedure :: read_markers => trdr_read_markers
end type trajectory_reader_t
 
contains

!******************************************************************************

subroutine trdr_init(this, fn_fns, col_nums, ierr, from_file)
    !!  Creates a `trajectory_reader_t` object.

    class(trajectory_reader_t), intent(out) :: this
        !! *trajectory_reader_t* instance.
    character(len=*), intent(in) :: fn_fns
        !! Names of the file containing names of dump files. Alternatively,
        !! if `from_file` is set to true, `fn_fns` contains all file names as
        !! well as the frame markers.
    integer, dimension(:), intent(in) :: col_nums
        !! Column numbers to be read from dump files
    integer, intent(out) :: ierr
        !! Error flag {0, 1}
    logical, intent(in), optional :: from_file
        !! Whether to read in frame markers from file `fn_fns`.
    character(len=:), allocatable :: fn
    integer(ip_long), dimension(:,:), allocatable :: tmp
    integer(ip_long), dimension(:), allocatable :: markers
    integer :: i, ios, fu, n, m, mxfnln, nframes, cursize, newsize 

    ierr = 0
    this%col_nums = col_nums

    !If markers are to be read from file
    if (present(from_file) .and. from_file) then
        call this%read_markers(fn_fns)
        return
    end if

    !Initial allocation, will be enlarged if necessary
    allocate( this%frame_markers(2,2) )
    cursize = 2 !Current size

    !Read in dump file names. First pass: get max file name length; 
    !second pass: read the file names
    open(newunit=fu, file=fn_fns, action='read', status='old')
    read(fu, *) n
    mxfnln = 0
    do i = 1, n
        call readline(fu, fn, '', ios)
        if (ios /= 0) exit
        mxfnln = max(mxfnln, len(fn))
    end do
    close(fu)

    allocate(character(len=mxfnln):: this%fns_dump(n))

    open(newunit=fu, file=fn_fns, action='read', status='old')
    read(fu, *) n
    do i = 1, n
        read(fu, "(a)") this%fns_dump(i)
    end do
    close(fu)

    !Get frame markers
    nframes = 0; cursize = 2
    do i = 1, size(this%fns_dump)
        fn = this%fns_dump(i)
        if (str_endswith(fn, '.bin')) then
            call binary_get_frame_markers(fn, m, markers) !m is the number of frames
        else
            call ascii_get_frame_markers(fn, m, markers)
        end if
        if (m == 0) cycle !No frames in this file
        !Enlarge this%frame_markers if necessary
        if (cursize < (nframes+m)) then
            newsize = 2*(nframes+m)
            allocate(tmp(2,newsize))
            tmp(:,1:cursize) = this%frame_markers
            call move_alloc(tmp, this%frame_markers)
            cursize = newsize
        end if
        this%frame_markers(1,nframes+1:nframes+m) = i
        this%frame_markers(2,nframes+1:nframes+m) = markers(1:m)
        nframes = nframes + m
        write(*,"(a,1x,a,1x,a)") "Dump file", fn, "marker set"
    end do

    !Shrink to fit
    if (cursize > nframes) then
        allocate(tmp(2,nframes))
        tmp = this%frame_markers(:,1:nframes)
        call move_alloc(tmp, this%frame_markers)
    end if

    this%num_frames = nframes

    end subroutine

!******************************************************************************

subroutine trdr_delete(this)
   !! After a call to this subroutine, all memory within `this` is deallocated,
   !! all components of `this` are reset to zero, and any underlying files are
   !! closed (if open).

    class(trajectory_reader_t), intent(in out) :: this

    this%num_frames = 0

    if (allocated(this%col_nums)) deallocate(this%col_nums)
    if (allocated(this%fns_dump)) deallocate(this%fns_dump)
    if (allocated(this%frame_markers)) deallocate(this%frame_markers)

    end subroutine

!******************************************************************************

subroutine trdr_write_markers(this, fn)
    !! Write dump file names and markers for each frame to a file.

    class(trajectory_reader_t), intent(in) :: this
    character(*), intent(in) :: fn
    integer :: fu, i, n

    open(newunit=fu, file=fn, action='write', status='replace')

    n = size(this%fns_dump)
    write(fu, "(a,1x,i0)") 'DUMP_FILES', n
    do i = 1, n
        write(fu, "(a)") this%fns_dump(i)
    end do

    n = size(this%frame_markers, 2)
    write(fu, "(a,1x,i0)") 'MARKERS', n
    do i = 1, n
        write(fu, "(i0,1x,i0)") this%frame_markers(:,i)
    end do

    close(fu)

    end subroutine

!******************************************************************************

subroutine trdr_read_markers(this, fn)
    !! Read dump file names and markers for each frame from a file.

    class(trajectory_reader_t), intent(in out) :: this
    character(*), intent(in) :: fn
    character(:), allocatable :: line, words(:)
    integer :: fu, i, ios, n, mxfnln

    open(newunit=fu, file=fn, action='read', status='old')
    do
        call readline(fu, line, '', ios)
        if (ios /= 0) exit
        if (str_startswith(line, 'DUMP_FILES')) then
            call str_split_to_array(line, words)
            read(words(2),*) n
            mxfnln = 0
            do i = 1, n
                call readline(fu, line, '', ios)
                if (ios /= 0) exit
                mxfnln = max(mxfnln, len(line))
            end do
            if (allocated(this%fns_dump)) deallocate(this%fns_dump)
            allocate(character(len=mxfnln):: this%fns_dump(n))
        end if
    end do
    close(fu)

    open(newunit=fu, file=fn, action='read', status='old')
    do
        call readline(fu, line, '', ios)
        if (ios /= 0) exit
        if (str_startswith(line, 'DUMP_FILES')) then
            n = size(this%fns_dump)
            do i = 1, n
                read(fu, "(a)") this%fns_dump(i)
            end do
        end if

        if (str_startswith(line, 'MARKERS')) then
            call str_split_to_array(line, words)
            read(words(2),*) n
            if (allocated(this%frame_markers)) deallocate(this%frame_markers)
            allocate(this%frame_markers(2,n))
            do i = 1, n
                read(fu, *) this%frame_markers(:,i)
            end do
            this%num_frames = n
        end if
    end do
    close(fu)

    end subroutine

!******************************************************************************

subroutine trdr_read_frame(this, iframe, nts, natoms, is_triclinic, &
        boundary, box, dmp_data, ierr)
    !! Read a frame from an open trajectory.

    class(trajectory_reader_t), intent(in out) :: this
    integer, intent(in) :: iframe
        !! Frame number
    integer(ip_long), intent(out) :: nts
        !! Time step counter
    integer(ip_long), intent(out) :: natoms
        !! Number of atoms in the frame
    logical, intent(out) :: is_triclinic
        !! Is the simulation box triclinic? {T,F}
    character(len=2), dimension(3), intent(out) :: boundary
        !! Nature of the boundary (periodic, shrink wrap, etc...)
    real(rp), dimension(9), intent(out) :: box
        !! Box dimensions: xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz
    real(rp), dimension(:,:), allocatable, intent(in out) :: dmp_data
        !! (*m*,*n*) array, where `m >= len(col_nums)` and `n >= natoms`.
    integer, intent(out) :: ierr
        !! Error flag {0,1}.
    character(len=:), allocatable :: fn_dump
    integer(ip_long) :: frm_beg

    ierr = 0
    if ( (iframe < 1) .or. (iframe > this%num_frames) ) then
        write(*,"(a,i0,a,i0)") "ERROR: `iframe` out of bounds. &
            &`iframe` must be >= 1 and <= ", this%num_frames, &
            ". `iframe` = ", iframe
        ierr = 1; return
    end if

    fn_dump = this%fns_dump(this%frame_markers(1,iframe))
    frm_beg = this%frame_markers(2,iframe)

    if (str_endswith(fn_dump, '.bin')) then
        call binary_read_frame(fn_dump, frm_beg, this%col_nums, nts, natoms, &
            is_triclinic, boundary, box, dmp_data)
    else
        call ascii_read_frame(fn_dump, frm_beg, this%col_nums, nts, natoms, &
            is_triclinic, boundary, box, dmp_data)
    end if
        
    end subroutine

!******************************************************************************

subroutine ascii_get_frame_markers(fn_dump, num_frames, frame_markers)
    !! Determines the total number of frames and the starting line number of each
    !! frame in a LAMMPS ascii dump file.

    character(len=*), intent(in) :: fn_dump
        !! Name of the ascii dump file.
    integer, intent(out) :: num_frames
        !! Number of frames in `fn_dump`.
    integer(ip_long), dimension(:), allocatable, intent(in out) :: frame_markers
        !! (`num_frames`,) array of the starting line number of each frame in
        !! `fn_dump`.
    character(len=:), allocatable :: line
    integer(ip_long), dimension(:), allocatable :: tmp
    integer :: fu, nframes, natoms, lineno, lineno_beg, i, ios
    
    nframes = 0; lineno = 0
    if (allocated(frame_markers)) then
        frame_markers = 0
    else
        !Initial allocation, may be increased later if necessary
        allocate(frame_markers(2)) 
        frame_markers = 0
    end if

    open(newunit=fu, file=fn_dump, action='read', status='old')

    do 
        call readline(fu, line, '', ios)
        if (ios /= 0) exit
        lineno = lineno + 1

        if (str_startswith(line, 'ITEM: TIMESTEP')) then
            lineno_beg = lineno
            call readline(fu, line, '', ios)
            if (ios /= 0) exit
            lineno = lineno + 1
        end if
        if (str_startswith(line, 'ITEM: NUMBER OF ATOMS')) then
            read(fu, *, iostat=ios) natoms
            if (ios /= 0) exit
            lineno = lineno + 1
        end if
        if (str_startswith(line, 'ITEM: ATOMS')) then
            do i = 1, natoms
                call readline(fu, line, '', ios)
                if (ios /= 0) exit
                lineno = lineno + 1
            end do

            nframes = nframes + 1
            if (size(frame_markers) < nframes) then
                if (allocated(tmp)) deallocate(tmp)
                allocate(tmp(2*size(frame_markers)))
                tmp(:size(frame_markers)) = frame_markers
                call move_alloc(tmp, frame_markers)
            end if
            frame_markers(nframes) = lineno_beg
        end if
    end do

    num_frames = nframes
    close(fu)

    end subroutine

!******************************************************************************

subroutine ascii_read_frame(fn_dump, frm_beg, col_nums, nts, natoms, &
    is_triclinic, boundary, box, dmp_data)
    !! Reads a frame from a LAMMPS ascii dump file.

    character(len=*), intent(in) :: fn_dump
        !! Name of the ascii dump file
    integer(ip_long), intent(in) :: frm_beg
        !! Line number corresponding to the beginning of the frame
    integer, dimension(:), intent(in) :: col_nums
        !! Column numbers corresponding to data to be retrieved
    integer(ip_long), intent(out) :: nts
        !! Time step counter
    integer(ip_long), intent(out) :: natoms
        !! Number of atoms in the frame
    logical, intent(out) :: is_triclinic
        !! Is the simulation box triclinic? {T,F}
    character(len=2), dimension(3), intent(out) :: boundary
        !! Nature of the boundary (periodic, shrink wrap, etc...)
    real(rp), dimension(9), intent(out) :: box
        !! Entries are: xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz. If 
        !! `is_triclinic`  is false, xy, xz, and yz will be set to zero.
    real(rp), dimension(:,:), allocatable, intent(in out) :: dmp_data
        !! (*m*,*n*) array, where `m >= len(col_nums)` and `n >= natoms`.
    character(len=:), dimension(:), allocatable :: words
    character(len=:), allocatable :: line
    integer :: fu, i, j, ios, nrows, ncols
    logical :: frame_is_read
    
    frame_is_read = .false.
    box = 0.0_rp

    open(newunit=fu, file=fn_dump, action='read', status='old')
    !Skip to the beginning of the frame
    do i = 1, int(frm_beg-1, ip)
        read(fu,*)
    end do

    do 
        call readline(fu, line, '', ios)
        if (frame_is_read) exit

        if (str_startswith(line, 'ITEM: TIMESTEP')) read(fu,*) nts

        if (str_startswith(line, 'ITEM: NUMBER OF ATOMS')) then
            read(fu,*) natoms
            if (allocated(dmp_data)) then
                nrows = size(dmp_data,1); ncols = size(dmp_data,2)
                if ( (nrows < size(col_nums)) .or. (ncols < natoms) ) then
                    deallocate(dmp_data)
                    allocate( dmp_data(size(col_nums), natoms) )
                end if
            else
                allocate( dmp_data(size(col_nums), natoms) )
            end if
            dmp_data = 0.0_rp
        end if

        if (str_startswith(line, 'ITEM: BOX BOUNDS')) then
            call str_split_to_array(line, words)
            boundary(1) = trim(words(4))
            boundary(2) = trim(words(5))
            boundary(3) = trim(words(6))
            if (size(words)==6) then
                read(fu,*) box(1), box(2)
                read(fu,*) box(3), box(4)
                read(fu,*) box(5), box(6)
                is_triclinic = .false.
            else
                read(fu,*) box(1), box(2), box(7)
                read(fu,*) box(3), box(4), box(8)
                read(fu,*) box(5), box(6), box(9)
                is_triclinic = .true.
            end if
        end if

        if (str_startswith(line, 'ITEM: ATOMS')) then
            call str_split_to_array(line, words)
            do i = 1, int(natoms, ip)
                call readline(fu, line, '', ios)
                call str_split_to_array(line, words)
                do j = 1, size(col_nums)
                    read(words(col_nums(j)),*) dmp_data(j, i)
                end do
            end do
            frame_is_read = .true.
        end if
    end do

    close(fu)

    end subroutine

!******************************************************************************

subroutine binary_get_frame_markers(fn_dump, num_frames, frame_markers)
    !! Determines the total number of frames and the starting position of each
    !! frame from a LAMMPS binary dump file.

    character(len=*), intent(in) :: fn_dump
        !! Name of the binary dump file.
    integer, intent(out) :: num_frames
        !! Number of frames in `fn_dump`.
    integer(ip_long), dimension(:), allocatable, intent(out) :: frame_markers
        !! (`num_frames`,) array of the starting position of each frame in
        !! `fn_dump`.
    integer(ip_long), dimension(:), allocatable :: tmp

    integer(ip_long) :: ntimestep, natoms, len_ms
    integer :: nchunk, triclinic, len_us, len_cols, len_buf, ichunk, size_one
    integer :: endian, revision, revold
    real(rp) :: xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz, time
    integer, dimension(2,3) :: boundary
    real(rp), dimension(:), allocatable :: buffer
    character(len=:), allocatable :: magic_string, columns, unit_style
    character(len=1) :: flag
    integer(ip_long) :: pos, pos_beg
    integer :: fu, nframes, ios
    data       endian, revision, revold /Z'0001', Z'0001', Z'0001'/

    
    nframes = 0; pos = 0
    if (allocated(frame_markers)) then
        frame_markers = 0
    else
        !Initial allocation, may be increased later if necessary
        allocate(frame_markers(2)) 
        frame_markers = 0
    end if

    open(newunit=fu, file=fn_dump, access='stream', form='unformatted', &
        action='read', status='old')

    do 
        read(fu, iostat=ios) ntimestep
        if (ios /= 0) exit
        pos_beg = pos + 1
        pos = pos + sizeof_long_int

        if (ntimestep < 0) then
            !First int64 encodes negative of the length of the magic string
            len_ms = -ntimestep
            if (allocated(magic_string)) deallocate(magic_string)
            allocate(character(len=len_ms) :: magic_string)

            read(fu, iostat=ios) magic_string
            if (ios /= 0) exit
            pos = pos + len_ms*sizeof_char

            !Read endian flag
            read(fu, iostat=ios) endian
            if (ios /= 0) exit
            pos = pos + sizeof_int

            !Read revision number
            read(fu, iostat=ios) revision
            if (ios /= 0) exit
            pos = pos + sizeof_int

            !Read actual value of ntimestep
            read(fu, iostat=ios) ntimestep
            if (ios /= 0) exit
            pos = pos + sizeof_long_int
        end if

        read(fu, iostat=ios) natoms, triclinic, boundary, xlo, xhi, &
                            ylo, yhi, zlo, zhi
        if (ios /= 0) exit
        pos = pos + sizeof_long_int + sizeof_int + 6*sizeof_int + 6*sizeof_real
        if (triclinic /= 0) then
            read(fu, iostat=ios) xy, xz, yz
            if (ios /= 0) exit
            pos = pos + 3*sizeof_real
        end if

        read(fu, iostat=ios) size_one 
        if (ios /= 0) exit
        pos = pos + sizeof_int

        if ( (len(magic_string) > 0) .and. (revision > revold) ) then
            !Newer format includes units string, columns string & time
            read(fu, iostat=ios) len_us
            if (ios /= 0) exit
            pos = pos + sizeof_int

            if (len_us > 0) then
                if (allocated(unit_style)) deallocate(unit_style)
                allocate(character(len=len_us):: unit_style)
                read(fu, iostat=ios) unit_style
                if (ios /= 0) exit
                pos = pos + len_us*sizeof_char
            end if

            read(fu, iostat=ios) flag
            if (ios /= 0) exit
            pos = pos + sizeof_char
            if (iachar(flag) /= 0) then
                read(fu, iostat=ios) time
                if (ios /= 0) exit
                pos = pos + sizeof_real
            end if

            read(fu, iostat=ios) len_cols
            if (ios /= 0) exit
            pos = pos + sizeof_int

            if (allocated(columns)) deallocate(columns)
            allocate(character(len=len_cols):: columns)
            read(fu, iostat=ios) columns
            if (ios /= 0) exit
            pos = pos + len_cols*sizeof_char
        end if

        read(fu, iostat=ios) nchunk
        if (ios /= 0) exit
        pos = pos + sizeof_int
        !Loop over processor chunks
        do ichunk = 1, nchunk
            read(fu, iostat=ios) len_buf
            if (ios /= 0) exit
            pos = pos + sizeof_int
            if (.not. allocated(buffer)) allocate(buffer(len_buf))
            if (size(buffer) < len_buf) then
                deallocate(buffer)
                allocate(buffer(len_buf))
            end if
            !Read chunk
            read(fu, iostat=ios) buffer(1:len_buf)
            if (ios /= 0) exit
            pos = pos + len_buf*sizeof_real
        end do

        nframes = nframes + 1
        if (size(frame_markers) < nframes) then
            if (allocated(tmp)) deallocate(tmp)
            allocate(tmp(2*size(frame_markers)))
            tmp(:size(frame_markers)) = frame_markers
            call move_alloc(tmp, frame_markers)
        end if
        frame_markers(nframes) = pos_beg
    end do

    num_frames = nframes
    close(fu)

    end subroutine

!******************************************************************************

subroutine binary_read_frame(fn_dump, frm_beg, col_nums, nts, natoms, &
    is_triclinic, boundary, box, dmp_data)
    !! Reads a frame from a LAMMPS binary dump file.

    character(len=*), intent(in) :: fn_dump
        !! Name of the binary dump file
    integer(ip_long), intent(in) :: frm_beg
        !! Line number corresponding to the beginning of the frame
    integer, dimension(:), intent(in) :: col_nums
        !! Column numbers corresponding to data to be retrieved
    integer(ip_long), intent(out) :: nts
        !! Time step counter
    integer(ip_long), intent(out) :: natoms
        !! Number of atoms in the frame
    logical, intent(out) :: is_triclinic
        !! Is the simulation box triclinic? {T,F}
    character(len=2), dimension(3), intent(out) :: boundary
        !! Nature of the boundary (periodic, shrink wrap, etc...)
    real(rp), dimension(9), intent(out) :: box
        !! Entries are: xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz. If 
        !! `is_triclinic`  is false, xy, xz, and yz will be set to zero.
    real(rp), dimension(:,:), allocatable, intent(in out) :: dmp_data
        !! (*m*,*n*) array, where `m >= len(col_nums)` and `n >= natoms`.
    integer(ip_long) :: ntimestep
        ! Time step for the current frame
    integer(ip_long) :: len_ms
        ! Length of the magic string (if present)
    integer :: len_us, len_cols, len_buf
        ! Lengths of the arrays unit_style, columns, and buffer.
    integer :: size_one
        ! The number of values (real) per line
    integer :: nchunk
        ! Number of processor chunks
    integer :: triclinic
        ! Flag for triclinic box {0, 1}
    real(rp) :: xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz
        ! Box dimensions, tilt factors, etc. for triclinic boxes
    real(rp) :: time
        ! Elapsed simulation time (if present)
    integer, dimension(2,3) :: bndry
        ! Flags for boundry conditions as read from a binary file
    character(len=1) :: flag
        ! Flag indicating whether elapsed time data is present in the file
    character(len=:), allocatable :: magic_string
        ! String for new file format
    character(len=:), allocatable :: columns
        ! Names of per-atom data fields
    character(len=:), allocatable :: unit_style
        ! LAMMPS unit field (if present)
    real(rp), dimension(:), allocatable :: buffer
        ! Per atom data buffer
    integer :: endian, revision, revold
    integer :: i, j, k, iline, ichunk, ibeg, ios, fu, nrows, ncols, nlines
    data       endian, revision, revold /Z'0001', Z'0001', Z'0001'/
    
    box = 0.0_rp

    open(newunit=fu, file=fn_dump, access='stream', form='unformatted', &
        action='read', status='old')

    read(fu, pos=frm_beg, iostat=ios) ntimestep
    !Detect whether new or old format
    if (ntimestep < 0) then
        !New format: First int64 encodes negative of the length of the magic string
        len_ms = -ntimestep
        if (allocated(magic_string)) deallocate(magic_string)
        allocate(character(len=len_ms)::magic_string)
        read(fu) magic_string
        !Read endian flag, revision number, and actual value of `ntimestep`
        read(fu) endian, revision, ntimestep
    else
        !Old format: No magic string, etc.
        magic_string = ''
    end if
    nts = ntimestep

    read(fu) natoms, triclinic, bndry, xlo, xhi, ylo, yhi, zlo, zhi
    box(1) = xlo; box(2) = xhi
    box(3) = ylo; box(4) = yhi
    box(5) = zlo; box(6) = zhi

    if (triclinic /= 0) then
        is_triclinic = .true.
        read(fu) xy, xz, yz
        box(7) = xy; box(8) = xz; box(9) = yz
    else
        is_triclinic = .false.
    end if

    do j = 1, 3
        do i = 1, 2
            select case ( int(bndry(i,j)) )
            case (0)
                boundary(j)(i:i) = 'p'
            case (1)
                boundary(j)(i:i) = 'f'
            case (2)
                boundary(j)(i:i) = 's'
            case (3)
                boundary(j)(i:i) = 'm'
            end select
        end do
    end do

    read(fu) size_one

    if (allocated(dmp_data)) then
        nrows = size(dmp_data,1); ncols = size(dmp_data,2)
        if ( (nrows < size(col_nums)) .or. (ncols < natoms) ) then
            deallocate(dmp_data)
            allocate( dmp_data(size(col_nums), natoms) )
        end if
    else
        allocate( dmp_data(size(col_nums), natoms) )
    end if
    dmp_data = 0.0_rp

    if ( (len(magic_string) > 0) .and. (revision > revold) ) then
        !Newer format includes units string, columns string & time
        read(fu) len_us
        if (len_us > 0) then
            if (allocated(unit_style)) deallocate(unit_style)
            allocate(character(len=len_us):: unit_style)
            read(fu) unit_style
        end if

        read(fu) flag
        if (iachar(flag) /= 0) then
            read(fu) time
        end if

        read(fu) len_cols
        if (allocated(columns)) deallocate(columns)
        allocate(character(len=len_cols):: columns)
        read(fu) columns
    end if

    !Loop over processor chunks
    read(fu) nchunk
    k = 0
    do ichunk = 1, nchunk
        read(fu) len_buf
        if (.not. allocated(buffer)) allocate(buffer(len_buf))
        if (size(buffer) < len_buf) then
            deallocate(buffer)
            allocate(buffer(len_buf))
        end if
        !Read chunk
        read(fu) buffer(1:len_buf)
        nlines = len_buf/size_one

        do iline = 1, nlines
            ibeg = (iline-1)*size_one
            k = k + 1
            do j = 1, size(col_nums)
                dmp_data(j,k) = buffer(ibeg+col_nums(j))
            end do
        end do
    end do

    close(fu)

    end subroutine

!******************************************************************************

end module trajectory_reader_m
