!********************************************************************************!
!                                                                                !
! The MIT License (MIT)                                                          !
!                                                                                !
! Copyright (c) 2022 Sarit Dutta <saridut@gmail.com>                             !
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

module lcnpl_m
    !! Bookkeeping arrays for ligand coated NPLs.
 
use constants_m

implicit none

private
public :: lcnpl_t

type lcnpl_t
    integer(ip_long) :: num_atoms = 0

    contains
        procedure :: init       => trdr_init
        procedure :: delete     => trdr_delete
        procedure :: read_frame => trdr_read_frame
end type lcnpl_t
 
contains

!******************************************************************************

subroutine trdr_init(this, fns, ierr)
    !!  Creates a `trajectory_reader_t` object.

    class(trajectory_reader_t), intent(out) :: this
        !! *trajectory_reader_t* instance.
    character(len=*), dimension(:), intent(in) :: fns
        !! Names of the trajectory files.
    integer, intent(out) :: ierr
        !! Error flag {0, 1}
    character(len=:), allocatable :: fn_traj
    integer(ip_long) :: buf_pad_size
    integer :: ntraj, i
    
    ierr = 0
    this%fns = fns
    ntraj = size(fns)
    allocate( this%trajs(ntraj) )

    this%num_frames = huge(0); buf_pad_size = 0_ip_long
    do i = 1, ntraj
        !Open file for reading trajectory
        fn_traj = trim(fns(i))
        write(*, '(a,i0," : ",a)') "Traj ", i, fn_traj
        call this%trajs(i)%open(fn_traj, ierr)
        if (ierr /=0 ) return
        this%num_frames = min(this%num_frames, this%trajs(i)%num_frames)
        this%num_atoms = this%num_atoms + this%trajs(i)%num_atoms
        buf_pad_size = buf_pad_size &
                        + int(this%trajs(i)%buf_pad_size/sizeof_real, ip_long)
    end do
    this%num_fields = this%trajs(1)%size_one !from the first file, all same
    allocate(this%buf_pad(buf_pad_size))

    end subroutine

!******************************************************************************

subroutine trdr_delete(this)
   !! After a call to this subroutine, all memory within `this` is deallocated,
   !! all components of `this` are reset to zero, and any underlying files are
   !! closed (if open).

    class(trajectory_reader_t), intent(in out) :: this
    integer :: i

    do i = 1, size(this%trajs)
        call this%trajs(i)%close()
    end do

    if (allocated(this%trajs)) deallocate(this%trajs)
    if (allocated(this%fns)) deallocate(this%fns)
    if (allocated(this%buf_pad)) deallocate(this%buf_pad)

    this%num_fields = 0
    this%num_atoms = 0
    this%num_frames = 0

    end subroutine

!******************************************************************************

subroutine trdr_read_frame(this, iframe, nts, box, boundary, ierr)
    !! Read a frame from an open trajectory.

    class(trajectory_reader_t), intent(in out) :: this
    integer, intent(in)  :: iframe
        !! Frame number
    integer(ip_long), intent(out) :: nts
        !! Time step counter
    real(rp), dimension(9), intent(out) :: box
        !! Box dimensions: xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz
    integer, dimension(2,3), intent(out) :: boundary
        !! Boundary conditions
    integer, intent(out) :: ierr
        !! Error flag
    integer(ip_long) :: ibeg, iend
    integer :: i

    ierr = 0; ibeg = 1

    do i = 1, size(this%trajs)
        iend = ibeg + int(this%trajs(i)%buf_pad_size/sizeof_real, ip_long) - 1
        call this%trajs(i)%read(iframe, nts, box, boundary, &
            this%buf_pad(ibeg:iend), ierr)
        ibeg = iend + 1
    end do
        
    end subroutine

!******************************************************************************

end module lcnpl_m
