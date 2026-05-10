!********************************************************************************!
! The MIT License (MIT)                                                          !
!                                                                                !
! Copyright (c) 2023 Sarit Dutta <saridut@gmail.com>                             !
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
!********************************************************************************!

module simbox_m
!! Implements a simulation box with appropriate boundary conditions.

use constants_m

implicit none

type smbx_t
    character(len=2), dimension(3) :: boundary
    real(rp), dimension(3,3) :: basis
    real(rp), dimension(3,3) :: dl_basis
    real(rp) :: xlo, xhi, ylo, yhi, zlo, zhi
    real(rp) :: xy, xz, yz
    real(rp) :: volume
    logical :: is_deforming
    logical :: is_aligned
    logical :: is_triclinic
    contains
        procedure :: init => smbx_init
        procedure :: set_boundary => smbx_set_boundary
        procedure :: set_bounds => smbx_set_bounds
        procedure :: set_tilt => smbx_set_tilt
        procedure :: update => smbx_update
        procedure :: freeze => smbx_freeze
        procedure :: unfreeze => smbx_unfreeze
        procedure :: get_image => smbx_get_image
        procedure :: wrap_all => smbx_wrap_all
        procedure :: scale_all => smbx_scale_all
        procedure :: unscale_all => smbx_unscale_all
end type smbx_t

contains

!********************************************************************************

subroutine smbx_init(this, boundary, is_triclinic)
    !! Initializes an instance of *smbx_t*. Can also be called to reset.

    class(smbx_t), intent(out) :: this
        !! An instance of `smbx_t`.
    character, dimension(3), intent(in), optional :: boundary
    logical, intent(in), optional :: is_triclinic

    this%boundary = 'pp'; 
    this%xlo = 0.0_rp; this%ylo = 0.0_rp; this%zlo = 0.0_rp
    this%xhi = 0.0_rp; this%yhi = 0.0_rp; this%zhi = 0.0_rp
    this%xy = 0.0_rp; this%xz = 0.0_rp; this%yz = 0.0_rp
    this%basis = 0.0_rp; this%dl_basis = 0.0_rp; this%volume = 0.0_rp

    this%is_deforming = .false.
    this%is_aligned = .true.
    this%is_triclinic = .false.

    if (present(boundary)) this%boundary = boundary
    if (present(is_triclinic)) this%is_triclinic = is_triclinic

    end subroutine

!********************************************************************************

subroutine smbx_set_boundary(this, xb, yb, zb)
    !!Sets the boundary style of the box along the x, y, or z direction.

    class(smbx_t), intent(in out) :: this
    character(len=2), intent(in), optional :: xb
    character(len=2), intent(in), optional :: yb
    character(len=2), intent(in), optional :: zb

    if ( present(xb) ) this%boundary(1) = xb
    if ( present(yb) ) this%boundary(2) = yb
    if ( present(zb) ) this%boundary(3) = zb

    end subroutine

!********************************************************************************

subroutine smbx_set_bounds(this, xlo, xhi, ylo, yhi, zlo, zhi)
    !!Sets the upper or lower bound of the box along the x, y, or z direction.

    class(smbx_t), intent(in out)  :: this
    real(rp), intent(in), optional :: xlo
    real(rp), intent(in), optional :: xhi
    real(rp), intent(in), optional :: ylo
    real(rp), intent(in), optional :: yhi
    real(rp), intent(in), optional :: zlo
    real(rp), intent(in), optional :: zhi

    if ( present(xlo) ) this%xlo = xlo
    if ( present(xhi) ) this%xhi = xhi
    if ( present(ylo) ) this%ylo = ylo
    if ( present(yhi) ) this%yhi = yhi
    if ( present(zlo) ) this%zlo = zlo
    if ( present(zhi) ) this%zhi = zhi

    end subroutine

!********************************************************************************

subroutine smbx_set_tilt(this, xy, xz, yz)
    !!Sets the tilt for xy, xz, or yz.

    class(smbx_t), intent(in out) :: this
    real(rp), intent(in), optional :: xy
    real(rp), intent(in), optional :: xz
    real(rp), intent(in), optional :: yz

    if (present(xy)) this%xy = xy
    if (present(xz)) this%xz = xz
    if (present(yz)) this%yz = yz

    end subroutine

!********************************************************************************

subroutine smbx_update(this)
    !!Updates the basis vectors and volume.

    class(smbx_t), intent(in out) :: this
    real(rp), dimension(3) :: a, b, c

    this%basis = 0.0_rp
    this%basis(1,1) = this%xhi - this%xlo
    this%basis(2,2) = this%yhi - this%ylo
    this%basis(3,3) = this%zhi - this%zlo

    a = this%basis(:,1); b = this%basis(:,2); c = this%basis(:,3)

    this%volume = a(1)*( b(2)*c(3)-c(2)*b(3) ) - a(2)*( b(1)*c(3)-c(1)*b(3) ) &
            + a(3)*( b(1)*c(2)-c(1)*b(2) )

    end subroutine

!********************************************************************************

subroutine smbx_freeze(this)
    !!Specifies *this* as non-deforming.

    class(smbx_t), intent(in out) :: this

    this%is_deforming = .false.

    end subroutine

!********************************************************************************

subroutine smbx_unfreeze(this)
    !!Specifies *this* as deforming.

    class(smbx_t), intent(in out) :: this

    this%is_deforming = .true.

    end subroutine

!********************************************************************************

subroutine smbx_get_image(this, r)
    !!Returns the minimum image of *r* under PBC.
    !https://scicomp.stackexchange.com/questions/20165/periodic-boundary-conditions-for-triclinic-box

    class(smbx_t), intent(in) :: this
    real(rp), dimension(3), intent(in out) :: r
    real(rp), dimension(3) :: rf
    
    if (this%is_aligned) then
        if (this%boundary(1)(1:1)=='p') then
            r(1) = r(1) - this%basis(1,1)*nint( r(1)/this%basis(1,1) )
        end if

        if (this%boundary(2)(1:1)=='p') then
            r(2) = r(2) - this%basis(2,2)*nint( r(2)/this%basis(2,2) )
        end if

        if (this%boundary(3)(1:1)=='p') then
            r(3) = r(3) - this%basis(3,3)*nint( r(3)/this%basis(3,3) )
        end if
    else
        !TODO: fix this for triclinic box
        rf = matmul(this%dl_basis, r)
        rf = rf - nint(rf)
        r = matmul(this%basis, rf)
    end if
    
    end subroutine

!********************************************************************************

subroutine smbx_wrap_all(this, coords)
    !!Wraps atom positions w.r.t. periodic boundary conditions.
    !https://scicomp.stackexchange.com/questions/20165/periodic-boundary-conditions-for-triclinic-box

    class(smbx_t), intent(in) :: this
    real(rp), dimension(:,:), intent(in out) :: coords
    real(rp), dimension(3) :: rf
    real(rp), dimension(3) :: llb_corner ! Lower left back corner of the box
    real(rp), dimension(3) :: diag
    integer :: n, i
    
    n = size(coords,2)
    llb_corner = [this%xlo, this%ylo, this%zlo]

    do i = 1, n
        coords(:,i) = coords(:,i) - llb_corner
    end do

    if (this%is_aligned) then
        diag = [this%basis(1,1), this%basis(2,2), this%basis(3,3)]

        if (this%boundary(1)(1:1)=='p') then
            coords(1,:) = coords(1,:) - diag(1)*floor( coords(1,i)/diag(1) )
        end if

        if (this%boundary(2)(1:1)=='p') then
            coords(2,:) = coords(2,:) - diag(2)*floor( coords(2,i)/diag(2) )
        end if

        if (this%boundary(3)(1:1)=='p') then
            coords(3,:) = coords(3,:) - diag(3)*floor( coords(3,i)/diag(3) )
        end if

        !do i = 1, n
        !    coords(:,i) = coords(:,i) - diag*floor( coords(:,i)/diag )
        !end do
    else
        !TODO: fix this for triclinic box
        do i = 1, n
            rf = matmul(this%dl_basis, coords(:,i))
            rf = rf - floor(rf)
            coords(:,i) = matmul(this%basis, rf)
        end do
    end if

    do i = 1, n
        coords(:,i) = coords(:,i) + llb_corner
    end do

    end subroutine

!********************************************************************************

subroutine smbx_scale_all(this, coords, coord_desc)
    !!Scales atom positions with box dimensions.

    class(smbx_t), intent(in) :: this
    real(rp), dimension(:,:), intent(in out) :: coords
    character(len=2), dimension(3), intent(in) :: coord_desc
    real(rp), dimension(3) :: llb_corner ! Lower left back corner of the box
    real(rp), dimension(3) :: rdiag
    
    llb_corner = [this%xlo, this%ylo, this%zlo]
    rdiag = [1.0_rp/this%basis(1,1), 1.0_rp/this%basis(2,2), 1.0_rp/this%basis(3,3)]

    if (coord_desc(1)(1:1)=='u') coords(1,:) = (coords(1,:) - llb_corner(1))*rdiag(1)
    if (coord_desc(2)(1:1)=='u') coords(2,:) = (coords(2,:) - llb_corner(2))*rdiag(2)
    if (coord_desc(3)(1:1)=='u') coords(3,:) = (coords(3,:) - llb_corner(3))*rdiag(3)

    end subroutine

!********************************************************************************

subroutine smbx_unscale_all(this, coords, coord_desc)
    !!Unscales atom positions with box dimensions.

    class(smbx_t), intent(in) :: this
    real(rp), dimension(:,:), intent(in out) :: coords
    character(len=2), dimension(3), intent(in) :: coord_desc
    real(rp), dimension(3) :: llb_corner ! Lower left back corner of the box
    real(rp), dimension(3) :: diag
    
    llb_corner = [this%xlo, this%ylo, this%zlo]
    diag = [this%basis(1,1), this%basis(2,2), this%basis(3,3)]

    if (coord_desc(1)(1:1)=='s') coords(1,:) = llb_corner(1) + coords(1,:)*diag(1)
    if (coord_desc(2)(1:1)=='s') coords(2,:) = llb_corner(2) + coords(2,:)*diag(2)
    if (coord_desc(3)(1:1)=='s') coords(3,:) = llb_corner(3) + coords(3,:)*diag(3)

    end subroutine

!********************************************************************************

end module simbox_m
