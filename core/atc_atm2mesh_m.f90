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

module atc_atm2mesh_m

use iso_c_binding
use mkl_spblas
use constants_m
use table_m

implicit none

type(matrix_descr)    :: spm_descr
type(sparse_matrix_t) :: spm
integer :: spm_nrows, spm_ncols

contains

!*******************************************************************************

subroutine make_atme_tab(nodes, atm_pos, rcut, meth, tab)
    !! Creates a table assigning atoms to mesh nodes

    real(rp), dimension(:,:), intent(in) :: nodes
    real(rp), dimension(:,:), intent(in) :: atm_pos
    real(rp), intent(in) :: rcut
    character(len=*), intent(in) :: meth
    type(itable_t), intent(out) :: tab
    real(rp) :: rcutsq, rij_sq
    real(rp), dimension(3) :: ri, rj, rij
    integer :: nnodes, na
    integer :: i, j, ierr

    nnodes = size(nodes,2); na = size(atm_pos,2)
    rcutsq = rcut*rcut

    call itbl_init(tab, nnodes, ierr)

    !Loop over all pairs to build list
    do i = 1, nnodes
        ri = nodes(:,i)
        do j = 1, na
            rj = atm_pos(:,j)
            rij = rj - ri
            if ( any(abs(rij) > rcut) ) cycle
            rij_sq = sum(rij**2)
            if (rij_sq < rcutsq) call tab%append(i,j)
        end do
    end do

    !Sort all rows
    call tab%sort()

    !Release additional memory
    call tab%shrink_to_fit()

    end subroutine

!*******************************************************************************

subroutine build_spread_mat(nodes, atm_pos, rcut, meth, sup, fnc)
    !! Construct the speading matrix

    real(rp), dimension(:,:), intent(in) :: nodes
    real(rp), dimension(:,:), intent(in) :: atm_pos
    real(rp), intent(in) :: rcut
    character(len=*), intent(in) :: meth
    character(len=*), intent(in) :: sup ! {'s': sphere, 'd': disc}
    character(len=*), intent(in) :: fnc ! {'uni': Uniform, 'qsp': Quartic spline}
    real(rp), dimension(:), allocatable :: val
    integer, dimension(:), allocatable :: col
    integer, dimension(:), allocatable :: rowindx
    type(itable_t) :: tab
    integer, dimension(:), pointer :: ptr
    real(rp), dimension(3) :: ri, rj, rij
    real(rp) :: rij_mag
    integer :: i, j, k1, n, info

    call make_atme_tab(nodes, atm_pos, rcut, meth, tab)

    n = tab%get_size()
    allocate(val(n), col(n), rowindx(tab%num_rows+1))

    k1 = 1
    do i = 1, tab%num_rows
        ri = nodes(:,i)
        call tab%get_row(i, ptr)
        n = size(ptr)
        rowindx(i) = k1
        do j = 1, n
            col(k1) = ptr(j)
            rj = atm_pos(:,ptr(j))
            rij = rj - ri
            rij_mag = norm2(rij)
            val(k1) = floc(rij_mag, rcut, sup, fnc) !value here
            k1 = k1 + 1
        end do
    end do
    rowindx(tab%num_rows+1) = k1

    spm_nrows = size(nodes,2); spm_ncols = size(atm_pos,2)
    info = mkl_sparse_d_create_csr(spm, SPARSE_INDEX_BASE_ONE, spm_nrows, &
            spm_ncols, rowindx, rowindx(2), col, val)
    
    if (info /= SPARSE_STATUS_SUCCESS) then
        print *, '  mkl_sparse_d_create_csr: ', info
        stop
    end if

    !Matrix descriptor
    spm_descr%type = SPARSE_MATRIX_TYPE_GENERAL

    !Matrix hints
    info = mkl_sparse_set_mv_hint(spm, sparse_operation_non_transpose, &
        spm_descr, 100)

    !Analyse sparse matrix
    info = mkl_sparse_optimize(spm)

    end subroutine

!*******************************************************************************

subroutine delete_spread_mat()
    !! Destroy the speading matrix

    integer :: info

    info = mkl_sparse_destroy(spm)

    if (info /= SPARSE_STATUS_SUCCESS) then
        print *, '  mkl_sparse_destroy: ', info
        stop
    end if

    end subroutine

!*******************************************************************************

subroutine atm_to_mesh(atm_data, msh_data)
    !! Spread atomic property to the mesh nodes

    real(rp), dimension(:), intent(in) :: atm_data
    real(rp), dimension(:), intent(out) :: msh_data
    real(rp) :: alpha, beta
    integer :: info

    msh_data = 0.0_rp

    alpha = 1.0_rp; beta = 0.0_rp

    info = mkl_sparse_d_mv(sparse_operation_non_transpose, alpha, spm, spm_descr, &
        atm_data, beta, msh_data)

    if (info /= SPARSE_STATUS_SUCCESS) then
        print *, '  mkl_sparse_d_mv: ', info
        stop
    end if

    end subroutine

!*******************************************************************************

pure function floc(r, rcut, sup, fnc) result(res)
    !! Localization function: Quartic spline or uniform on a sphere or disc

    real(rp), intent(in) :: r
    real(rp), intent(in) :: rcut
    character(len=1), intent(in) :: sup !{'s': sphere, 'd': disc}
    character(len=*), intent(in) :: fnc !{'uni': Uniform, 'qsp': Quartic spline}
    real(rp) :: res
    real(rp) :: rrcut, rrw, trm, rrcut_cub
    real(rp), parameter :: pf_s = 105.0_rp/(16.0_rp*math_pi)
    real(rp), parameter :: pf_d = 5.0_rp/math_pi
    real(rp) :: pf

    if (fnc == 'uni') then
        if (sup == 's') then
            res = 3.0_rp/(4.0_rp*math_pi*rcut**3)
        else if (sup == 'd') then
            res = 1.0_rp/(math_pi*rcut*rcut)
        end if
    else if (fnc=='qsp') then
        rrcut = 1.0_rp/rcut
        if (sup=='s') then
            pf = pf_s*rrcut*rrcut*rrcut
        else if (sup=='d') then
            pf = pf_d*rrcut*rrcut
        end if
        rrw = r*rrcut; trm = 1.0_rp - rrw
        res = pf*(1.0_rp + 3.0_rp*rrw)*trm*trm*trm
    end if

    end function

!*******************************************************************************

end module atc_atm2mesh_m
