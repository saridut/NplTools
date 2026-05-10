!********************************************************************************!
!                                                                                !
! The MIT License (MIT)                                                          !
!                                                                                !
! Copyright (c) 2022 Sarit Dutta                                                 !
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

module atmcfg_m

use constants_m
use specs_m

implicit none

!******************************************************************************

type atmcfg_t
    !Particle configuration: Atoms
    integer :: num_atom_types = 0
        !! Number of *atom_type*s
    type(atm_specs_t), dimension(:), allocatable :: atom_specs
        !! (*num_atom_types*,) array.
    integer, dimension(:), allocatable :: atom_pop
        !! (*num_atom_types*,) array. Population of atoms of each type.
    integer :: num_atoms = 0
        !!  Number of atoms
    integer :: num_molecules = 0
        !!  Number of molecules
    integer , dimension(:,:), allocatable:: atoms
        !! (2, *num_atoms*,) array. For atom *i*, its type *at = atoms(1,i)*,
        !! molecule *(atoms(2,i))* charge *charge(i)*, position *coordinates(:,i)*, etc.
    real(rp), dimension(:), allocatable :: mass
        !! (*num_atoms*,) array.
    real(rp), dimension(:), allocatable :: charge
        !! (*num_atoms*,) array
    character(len=2), dimension(3) :: coord_desc = ['uw', 'uw', 'uw']
        !! Coordinate descriptor along each dimension. Allowed values
        !! are: `uw`: unscaled wrapped, `sw`: scaled wrapped, `uu`: unscaled
        !! unwrapped, and `su`: scaled unwrapped.
    real(rp), dimension(:,:), allocatable :: coordinates
        !! (3, *num_atoms*) array
    integer, dimension(:,:), allocatable :: image_box
        !! (3, *num_atoms*) array
    real(rp), dimension(:,:), allocatable :: velocities
        !!  (3, *num_atoms*) array
    real(rp), dimension(:,:), allocatable :: forces
        !!  (3, *num_atoms*) array
    real(rp), dimension(:,:), allocatable :: stresses
        !!  (6, *num_atoms*) array. Symmetric per-atom stress, where for atom
        !! *i*, the elements of *stresses(:,i)* denote components in the
        !! following order: `xx`, `yy`, `zz`, `xy`, `xz`, 'yz`.

    !Particle configuration: VDW (pair) interactions
    integer :: num_vdw_types = 0
        !!  Number of *vdw_type*s
    integer, dimension(:,:), allocatable :: vdw_pairs
        !!  (2, *num_vdw_types*) array. Stores atom type of interacting pairs, such
        !! that at_i >= at_j.
    type(ia_specs_t), dimension(:), allocatable :: vdw_specs
         !!  (*num_vdw_types*,) array.
    
    !Particle configuration: Bonds
    integer :: num_bond_types = 0
        !!  Number of *bond_type*s
    type(ia_specs_t), dimension(:), allocatable :: bond_specs
         !!  (*num_bond_types*,) array
    integer :: num_bonds = 0
        !!  Total number of bonds.
    integer, dimension(:,:), allocatable :: bonds
        !! (3, *num_bonds*) array. Bond *i* is of type *bt = bonds(1,i)*,  directed from
        !! atom *bonds(2,i)* to *bonds(3,i)*. 

    !Particle configuration: Angles
    integer :: num_angle_types = 0
        !!  Number of *angle_type*s
    type(ia_specs_t), dimension(:), allocatable :: angle_specs
         !!  (*num_angle_types*,) array
    integer :: num_angles = 0
        !!  Number of angles
    integer, dimension(:,:), allocatable :: angles
        !! (4, *num_angles*) array. Angle *i* is of type *ant = angles(1,i)*, incident
        !! to atoms *angles(2,i)*, *angles(3,i)*, and *angles(4,i)*.
    
    !Particle configuration: Dihedrals
    integer :: num_dihedral_types = 0
        !!  Number of *dihedral_type*s
    type(ia_specs_t), dimension(:), allocatable :: dihedral_specs
         !!  (*num_dihedral_types*,) array
    integer :: num_dihedrals = 0
        !!  Number of dihedrals
    integer, dimension(:,:), allocatable :: dihedrals
        !! (5, *num_dihedrals*) array. Dihedral *i* is of type
        !! *dht = dihedrals(1,i)*, incident to atoms *dihedrals(2,i)*,
        !! *dihedrals(3,i)*, and *dihedrals(4,i)*, and *dihedrals(5,i)*.

    !Particle configuration: Impropers
    integer :: num_improper_types = 0
        !!  Number of *improper_type*s
    type(ia_specs_t), dimension(:), allocatable :: improper_specs
         !!  (*num_improper_types*,) array
    integer :: num_impropers = 0
        !!  Number of impropers
    integer, dimension(:,:), allocatable :: impropers
        !! (5, *num_impropers*) array. Improper *i* is of type
        !! *imt = impropers(1,i)*, incident to atoms *impropers(2,i)*,
        !! *impropers(3,i)*, and *impropers(4,i)*, and *impropers(5,i)*.
end type atmcfg_t

contains

!******************************************************************************

subroutine atmcfg_delete(this)
    !! Deallocates all memory acquired by a `configuration_t` object and resets
    !! all other components to zero. Exception: `num_coeffs` is not reset to
    !! zero.

    type(atmcfg_t), intent(in out) :: this

    this%num_atom_types = 0; this%num_atoms = 0; this%num_molecules = 0
    this%coord_desc = ['uw', 'uw', 'uw']
    if (allocated(this%atom_specs))  deallocate(this%atom_specs)
    if (allocated(this%atom_pop  ))  deallocate(this%atom_pop)
    if (allocated(this%atoms     ))  deallocate(this%atoms)
    if (allocated(this%mass      ))  deallocate(this%mass)
    if (allocated(this%charge     )) deallocate(this%charge)
    if (allocated(this%coordinates)) deallocate(this%coordinates)
    if (allocated(this%image_box  )) deallocate(this%image_box)
    if (allocated(this%velocities )) deallocate(this%velocities)
    if (allocated(this%forces     )) deallocate(this%forces)
    if (allocated(this%stresses   )) deallocate(this%stresses)

    this%num_vdw_types = 0
    if (allocated(this%vdw_pairs )) deallocate(this%vdw_pairs)
    if (allocated(this%vdw_specs )) deallocate(this%vdw_specs)

    this%num_bond_types = 0; this%num_bonds = 0
    if (allocated(this%bond_specs)) deallocate(this%bond_specs)
    if (allocated(this%bonds      )) deallocate(this%bonds)

    this%num_angle_types = 0; this%num_angles = 0
    if (allocated(this%angle_specs)) deallocate(this%angle_specs)
    if (allocated(this%angles      )) deallocate(this%angles)

    this%num_dihedral_types = 0; this%num_dihedrals = 0
    if (allocated(this%dihedral_specs)) deallocate(this%dihedral_specs)
    if (allocated(this%dihedrals      )) deallocate(this%dihedrals)

    this%num_improper_types = 0; this%num_impropers = 0
    if (allocated(this%improper_specs)) deallocate(this%improper_specs)
    if (allocated(this%impropers      )) deallocate(this%impropers)

    end subroutine

!******************************************************************************

end module atmcfg_m
