!********************************************************************************!
! The MIT License (MIT)                                                          !
!                                                                                !
! Copyright (c) 2020 Sarit Dutta <saridut@gmail.com>                             !
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

module cell_list_m
    !! Sorts atoms using a cell list.
    !!
    !! The algorithm to build the cell list partially follows the techniques in Watanabe et
    !! al, 2011, “Efficient Implementations of Molecular Dynamics Simulations for
    !! Lennard-Jones Systems,” Prog. Theor. Phys. 126, 203–235.
    !!
    !! The pairlist is not explicitly built, rather the cells are directly
    !! looped over during force calculation.

use constants_m
use vector_m

implicit none

private

public cell_list_t

integer, parameter :: nncm = 13
    !! The 3D indices of a cell along with its neighboring cells are
    !! defined in terms of the offset array *nc_offset*, whose indices
    !! run from *-nncm* to *nncm*.
integer, dimension(3,2*nncm), parameter :: nc_offset = reshape( [ &
     -1,-1,-1,    0,-1,-1,    1,-1,-1, &
     -1, 1,-1,    0, 1,-1,    1, 1,-1, &
     -1, 0,-1,    1, 0,-1,    0, 0,-1, &
      0,-1, 0,    1,-1, 0,   -1,-1, 0, &
     -1, 0, 0,                1, 0, 0, &
      1, 1, 0,   -1, 1, 0,    0, 1, 0, &
      0, 0, 1,   -1, 0, 1,    1, 0, 1, &
     -1,-1, 1,    0,-1, 1,    1,-1, 1, &
     -1, 1, 1,    0, 1, 1,    1, 1, 1    ], [ 3, 2*nncm ] )

type cell_list_t
    logical :: is_half = .true.
        !! Whether half cell list is desired
    logical, dimension(3) :: is_periodic = .true.
        !! Is the boundary periodic along *x*, *y*, or *z* directions?
    real(rp), dimension(2,3) :: bounds = 0.0_rp
        !! Lower & upper bounds of the box
    real(rp) :: cell_size_nom = 0.0_rp
        !! Nominal cell size
    real(rp),  dimension(3) :: cell_size = 0.0_rp
        !! Actual cell size
    real(rp) :: fncmax = 1.25_rp
        !! Factor by which storage is incremented for any of the storage arrays.
        !! A value of 1.25 appears fine.
    integer,  dimension(3)   :: nc_max = 0
        !! Maximum number of cells along *x*, *y*, & *z*.
    integer,  dimension(3)   :: nc = 0
        !! Number of cells along *x*, *y*, & *z*.
    integer  :: nct_max = 0
        !! Maximum total number of cells.
    integer  :: nct = 0
        !! Total of cells.
    integer, dimension(:), allocatable :: cells
        !! *(na,)* array. Listing atoms in each cell.
    integer, dimension(:), allocatable :: cells_pos
        !! *(0:nct_max,)* index array. Note: 0-based indexing.
    type(ivector_t) :: cell_nbrs
        !! Lists neighbor cells for each cell.
    integer, dimension(:), allocatable :: cell_nbrs_pos
        !! *(0:nct_max,)* index array. Note: 0-based indexing.
    integer, dimension(:), allocatable :: host_cells
        !! *(na,)* array. *host_cells(i)* stores the linear index of the cell
        !! containing atom *i*. *na_max* is the total number of atoms under
        !! consideration.
    integer, dimension(:), allocatable :: cell_pop
        !! *(0:nct_max-1,)* array storing population of each cell. 
        !! Note: 0-based indexing.
    contains
        procedure :: init => cl_init
        procedure :: delete => cl_delete
        procedure :: print => cl_print
        procedure :: build => cl_build
        procedure :: update_box => cl_update_box
        procedure :: build_cell_nbrs => cl_build_cell_nbrs
        procedure :: get_num_cells => cl_get_num_cells
        procedure :: get_contents => cl_get_contents
        procedure :: get_nbr_cells => cl_get_nbr_cells
end type cell_list_t

contains

!******************************************************************************

subroutine cl_init(this, na, cs, bounds, is_periodic, is_half, fncmax)
    !! Initializes a cell list.

    class(cell_list_t), intent(out) :: this
    integer, intent(in) :: na
        !! Number of atoms to be handled.
    real(rp), intent(in) :: cs
        !! Nominal cell size.
    real(rp), dimension(2,3), intent(in) :: bounds
        !! Lower & upper bounds of the box
    logical, dimension(3), intent(in) :: is_periodic
        !! Is the boundary periodic along *x*, *y*, or *z* directions?
    logical, intent(in) :: is_half
        !! Half or full cell list?
    real(rp), intent(in), optional :: fncmax
    
    this%cell_size_nom = cs
    this%bounds = bounds
    this%is_periodic = is_periodic
    this%is_half = is_half
    if (present(fncmax)) this%fncmax = fncmax

    allocate( this%cells(na) )
    allocate( this%host_cells(na) )

    call this%build_cell_nbrs()

    end subroutine

!******************************************************************************

subroutine cl_update_box(this, bounds, is_periodic)
    !! Updates the bounding box for the cell list.

    class(cell_list_t), intent(in out) :: this
    real(rp), dimension(2,3), intent(in) :: bounds
        !! Lower & upper bounds of the box
    logical, dimension(3), optional, intent(in) :: is_periodic
        !! Is the boundary periodic along *x*, *y*, or *z* directions?

    this%bounds = bounds
    if (present(is_periodic)) this%is_periodic = is_periodic

    call this%build_cell_nbrs()

    end subroutine

!******************************************************************************

subroutine cl_build_cell_nbrs(this)
    !! Makes a table of neighboring cells.

    class(cell_list_t), intent(in out) :: this
    real(rp) :: len_box
    integer :: i, j
    integer :: ncx, ncy, ncz
    integer :: icx, icy, icz, ic
    integer :: jcx, jcy, jcz, jc
    integer :: ibeg, iend

    ! Sets the cell size. The actual cell size may be slightly larger.
    do i = 1, 3
        len_box = this%bounds(2,i)-this%bounds(1,i)
        this%nc(i) = floor( len_box / this%cell_size_nom )
        this%cell_size(i) = len_box / this%nc(i)
    end do
    this%nct = product(this%nc)

!   if ( (this%nct == 0) .and. (this%nct_max == 0) ) then
    if ( this%nct_max == 0 ) then
        !Buiding the nbr cell table for the first time
        this%nc_max = int(this%nc*this%fncmax)
        this%nct_max = product(this%nc_max)

        if (this%is_half) then
            call ivector_init(this%cell_nbrs, nncm*this%nct_max)
        else
            call ivector_init(this%cell_nbrs, 2*nncm*this%nct_max)
        end if

        allocate( this%cells_pos(0:this%nct_max) )
        allocate( this%cell_nbrs_pos(0:this%nct_max) )
        allocate( this%cell_pop(0:this%nct_max-1) )

    else if (this%nct > this%nct_max) then
        this%nc_max = int(this%nc*this%fncmax)
        this%nct_max = product(this%nc_max)

        deallocate(this%cells_pos)
        deallocate(this%cell_nbrs_pos)
        deallocate(this%cell_pop)

        call this%cell_nbrs%clear()

        allocate( this%cells_pos(0:this%nct_max) )
        allocate( this%cell_nbrs_pos(0:this%nct_max) )
        allocate( this%cell_pop(0:this%nct_max-1) )
    end if

    ncx = this%nc(1); ncy = this%nc(2); ncz = this%nc(3)

    if (this%is_half) then
        ibeg = nncm+1; iend = 2*nncm
    else
        ibeg = 1; iend = 2*nncm
    end if

    this%cell_nbrs_pos(0) = 1
    do icz = 0, ncz-1
        do icy = 0, ncy-1
            do icx = 0, ncx-1
                ic = icz*ncx*ncy + icy*ncx + icx
                this%cell_nbrs_pos(ic+1) = this%cell_nbrs_pos(ic)

                do j = ibeg, iend
                    jcx = icx + nc_offset(1,j)
                    jcy = icy + nc_offset(2,j)
                    jcz = icz + nc_offset(3,j)

                    if (this%is_periodic(1)) jcx = modulo(jcx, ncx)
                    if (this%is_periodic(2)) jcy = modulo(jcy, ncy)
                    if (this%is_periodic(3)) jcz = modulo(jcz, ncz)

                    if ( (jcx>=0) .and. (jcx<=ncx-1) .and. &
                         (jcy>=0) .and. (jcy<=ncy-1) .and. &
                         (jcz>=0) .and. (jcz<=ncz-1) ) then
                        jc = jcz*ncx*ncy + jcy*ncx + jcx
                        call this%cell_nbrs%append(jc)
                        this%cell_nbrs_pos(ic+1) = this%cell_nbrs_pos(ic+1) + 1
                    end if
                end do

            end do
        end do
    end do

    call this%cell_nbrs%shrink_to_fit()

    end subroutine

!******************************************************************************

subroutine cl_delete(this)
    !! Deallocates all memory associated with this.

    class(cell_list_t), intent(in out) :: this

    if (allocated(this%cells))         deallocate(this%cells)
    if (allocated(this%cells_pos))     deallocate(this%cells_pos)

    call this%cell_nbrs%delete()

    if (allocated(this%cell_nbrs_pos)) deallocate(this%cell_nbrs_pos)
    if (allocated(this%host_cells))    deallocate(this%host_cells)
    if (allocated(this%cell_pop))      deallocate(this%cell_pop)

    this%is_half = .true.; this%is_periodic = .true.; this%bounds = 0.0_rp
    this%cell_size_nom = 0.0_rp; this%cell_size = 0.0_rp; this%fncmax = 1.25_rp
    this%nc_max = 0; this%nc = 0
    this%nct_max = 0; this%nct = 0

    end subroutine

!******************************************************************************

subroutine cl_build(this, coords, atm_ids, ierr)
    !! Sorts atoms into cells

    class(cell_list_t), intent(in out) :: this
    real(rp), dimension(:,:), intent(in) :: coords
    integer, dimension(:), intent(in) :: atm_ids
    integer, intent(out) :: ierr
    real(rp), dimension(3) :: ri
    integer :: na, iatm, i, j
    integer :: ncx, ncy, ncz, nct
    integer :: icx, icy, icz, ic

    ierr = 0; na = size(atm_ids)

    !Check if storage needs to be increased
    if (na > size(this%cells)) then
        deallocate(this%cells); allocate(this%cells(na))
        deallocate(this%host_cells); allocate(this%host_cells(na))
    end if

    this%cells = 0; this%host_cells = 0; this%cell_pop = 0
    ncx = this%nc(1); ncy = this%nc(2); ncz = this%nc(3); nct = this%nct

    !Loop over atoms and put into cells
    do i = 1, na
        iatm = atm_ids(i)
        ri = coords(:,iatm)
        icx = int( ri(1)/this%cell_size(1) )
        icy = int( ri(2)/this%cell_size(2) )
        icz = int( ri(3)/this%cell_size(3) )

        if ((icx<0) .or. (icy<0) .or. (icz<0)) then
            write(*,'(a,i0)') 'Atom outside bounds. Id = ', iatm
            write(*,'(a,3(1x,es23.15))') 'Coordinates =', ri
            ierr = 1; return
        end if

        !If atoms are exactly on the box edge
        if ( icx > (ncx-1) ) icx = ncx - 1
        if ( icy > (ncy-1) ) icy = ncy - 1
        if ( icz > (ncz-1) ) icz = ncz - 1

        ic = icz*ncx*ncy + icy*ncx + icx
        this%host_cells(i) = ic
        this%cell_pop(ic) = this%cell_pop(ic) + 1
    end do

    !Loop over all cells to set up pointers
    this%cells_pos(0) = 1
    do ic = 0, (nct-1)
        this%cells_pos(ic+1) = this%cells_pos(ic) + this%cell_pop(ic)
    end do

    !Loop over all atoms
    do i = 1, na
        iatm = atm_ids(i)
        ic = this%host_cells(i)
        j = this%cells_pos(ic)
        this%cells(j) = iatm
        this%cells_pos(ic) = this%cells_pos(ic) + 1
    end do

    !Loop over all cells to set up pointers
    this%cells_pos(0) = 1
    do ic = 0, (nct-1)
        this%cells_pos(ic+1) = this%cells_pos(ic) + this%cell_pop(ic)
    end do

    end subroutine

!******************************************************************************

function cl_get_num_cells(this) result (res)
    !! Returns the total number of cells

    class(cell_list_t), intent(in) :: this
    integer :: res

    res = this%nct

    end function

!******************************************************************************

subroutine cl_get_contents(this, ic, res)
    !! Returns a pointer to the entries of cell with linear index *ic*.

    class(cell_list_t), intent(in), target :: this
    integer, intent(in) :: ic
    integer, dimension(:), pointer, intent(out) :: res
    integer :: ibeg, iend

    res => null()
    ibeg = this%cells_pos(ic); iend = this%cells_pos(ic+1) - 1
    res => this%cells(ibeg:iend)

    end subroutine

!******************************************************************************

subroutine cl_get_nbr_cells(this, ic, res)
    !! Returns a pointer to the neighbor cells of cell with linear index *ic*.

    class(cell_list_t), intent(in), target :: this
    integer, intent(in) :: ic
    integer, dimension(:), pointer, intent(out) :: res
    integer :: ibeg, iend

    res => null()
    ibeg = this%cell_nbrs_pos(ic); iend = this%cell_nbrs_pos(ic+1) - 1
    call this%cell_nbrs%get_data(res, ibeg, iend)

    end subroutine

!********************************************************************************

subroutine cl_print(this)
    !! Prints a cell list

    class(cell_list_t), intent(in), target :: this
    integer, dimension(:), pointer :: aic => null()
    integer, dimension(:), pointer :: nbrc => null()
    integer :: ncx, ncy, ncz
    integer :: icx, icy, icz, ic

    ncx = this%nc(1); ncy = this%nc(2); ncz = this%nc(3)

    write(*,'("ncx: ", i0, " ncy: ", i0, " ncz: ", i0)') ncx, ncy, ncz
    write(*,'("lcx: ", g0.6, " lcy: ", g0.6, " lcz: ", g0.6)') this%cell_size
    write(*, *) 'CELL CONTENTS'
    do icz = 0, ncz-1
        do icy = 0, ncy-1
            do icx = 0, ncx-1
                ic = icz*ncx*ncy + icy*ncx + icx
                call this%get_contents(ic, aic)
                if ( size(aic) > 0 ) then
                    write(*,'("(",i0,",",i0,",",i0,")[",i0,"] ")', advance='no') &
                        icx, icy, icz, ic
                    write(*, *) size(aic), aic
                end if
            end do
        end do
    end do

    write(*, *)
    write(*, *) 'NBR CELLS'
    do icz = 0, ncz-1
        do icy = 0, ncy-1
            do icx = 0, ncx-1
                ic = icz*ncx*ncy + icy*ncx + icx
                call this%get_nbr_cells(ic, nbrc)
                write(*,'("(",i0,",",i0,",",i0,")[",i0,"] ")', advance='no') &
                    icx, icy, icz, ic
                write(*, *) size(nbrc), nbrc
            end do
        end do
    end do

    end subroutine

!********************************************************************************

end module cell_list_m
