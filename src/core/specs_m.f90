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

module specs_m
!! **Type definitions for atoms, interactions, etc.**

use constants_m

implicit none

!******************************************************************************

type atm_specs_t
    character(len=8) :: name
    real(rp) :: mass
    integer  :: style !Style 0: Point particles, Style 1: Finite size particles.
end type atm_specs_t

!******************************************************************************

type ia_specs_t
    character(len=8) :: style
    real(rp), dimension(:), allocatable :: params
    logical, dimension(:), allocatable :: param_isint
        !! True if an parameter is an integer
    !! params(1) = vtltmin, params(2) = ftltmin
    !! params(3) = vtgtmax, params(4) = ftgtmax
    !! params(5) = intrp_meth. Interpolation method. 0: Linear, 1: Cubic spline
    !! params(6) = flag_bc. Boundary condition for cubic splines.
    !! 0: Natural (second derivative equals zero at the two ends),
    !! 1: Spline is quadratic over the first and last interval, 2: Not-a-knot condition
    real(rp), dimension(:), allocatable :: cache
    !! cache(1) = tmin, cache(2) = tmax
    integer  :: tab_size = 0
    real(rp), dimension(:), allocatable :: tab_t
        !! (*tab_size*,) array. Independent variable.
    real(rp), dimension(:), allocatable :: tab_v
        !! (*tab_size*,) array. Tabulated values of the potential.
    real(rp), dimension(:), allocatable :: tab_f
        !! (*tab_size*,) array. Tabulated values of the force.  Used only for
        !! linear interpolation.
    real(rp), dimension(:), allocatable :: tab_vpp
        !! (*tab_size*,) array. Stores tabulated second derivative of the potential.
        !! Used only for cubic spline interpolation.
end type ia_specs_t

!******************************************************************************

type extrn_field_t
    character(len=8) :: style
    real(rp), dimension(:), allocatable :: params
    real(rp), dimension(:), allocatable :: cache
end type extrn_field_t


!******************************************************************************

end module specs_m
