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

module config_io_m
    !! Routines for IO of config and dump files.

use constants_m
use strings_m
use simbox_m
use specs_m
use atmcfg_m

implicit none

contains

!******************************************************************************

subroutine read_ldf(simbox, atc, fn)
    !! Read from a LAMMPS data file

    type(smbx_t),     intent(out) :: simbox
    type(atmcfg_t),   intent(in out) :: atc
    character(len=*), intent(in)  :: fn
    integer :: atm_toff, bnd_toff, ang_toff, dhd_toff, imp_toff
    integer :: atm_idoff, bnd_idoff, ang_idoff, dhd_idoff, imp_idoff, &
                mol_idoff
    integer, dimension(:), allocatable    :: itmp
    integer, dimension(:,:), allocatable  :: itmp2
    real(rp), dimension(:), allocatable   :: rtmp
    real(rp), dimension(:,:), allocatable :: rtmp2
    type(atm_specs_t), dimension(:), allocatable :: asptmp
    type(ia_specs_t), dimension(:), allocatable :: isptmp
    real(rp) :: xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz
    character(:), allocatable :: line
    integer :: i, it, n, nt, at_i, imol, typ
    integer :: fu, ios
    logical :: is_triclinic = .false.

    atm_toff = atc%num_atom_types; atm_idoff = atc%num_atoms
    bnd_toff = atc%num_bond_types; bnd_idoff = atc%num_bonds
    ang_toff = atc%num_angle_types; ang_idoff = atc%num_angles
    dhd_toff = atc%num_dihedral_types; dhd_idoff = atc%num_dihedrals
    imp_toff = atc%num_improper_types; imp_idoff = atc%num_impropers
    mol_idoff = atc%num_molecules

    open(newunit=fu, file=fn, action='read', status='old')

    do 
        call readline(fu, line, '#', ios)
        if (ios /= 0) exit

        if (str_endswith(line, 'atom types')) then
            read(line,*) nt
            if (.not. allocated(atc%atom_specs))  then
                allocate(atc%atom_specs(nt))
                do it = 1, nt
                    atc%atom_specs(it)%name = str_from_num(it)
                    atc%atom_specs(it)%style = 0
                end do
            else
                allocate(asptmp(atm_toff+nt))
                asptmp(1:atm_toff) = atc%atom_specs
                do it = 1, nt
                    asptmp(atm_toff+it)%name = str_from_num(atm_toff+it)
                    asptmp(atm_toff+it)%style = 0
                end do
                call move_alloc(asptmp, atc%atom_specs)
            end if

            if (.not. allocated(atc%atom_pop))  then
                allocate(atc%atom_pop(nt))
                atc%atom_pop(:) = 0
            else
                allocate(itmp(atm_toff+nt))
                itmp = 0; itmp(1:atm_toff) = atc%atom_pop
                call move_alloc(itmp, atc%atom_pop)
            end if

            atc%num_atom_types = atm_toff + nt
        end if

        if (str_endswith(line, 'bond types')) then
            read(line,*) nt
            if (.not. allocated(atc%bond_specs))  then
                allocate(atc%bond_specs(nt))
            else
                allocate(isptmp(bnd_toff+nt))
                isptmp(1:bnd_toff) = atc%bond_specs
                call move_alloc(isptmp, atc%bond_specs)
            end if

            atc%num_bond_types = bnd_toff + nt
        end if

        if (str_endswith(line, 'angle types')) then
            read(line,*) nt
            if (.not. allocated(atc%angle_specs))  then
                allocate(atc%angle_specs(nt))
            else
                allocate(isptmp(ang_toff+nt))
                isptmp(1:ang_toff) = atc%angle_specs
                call move_alloc(isptmp, atc%angle_specs)
            end if

            atc%num_angle_types = ang_toff + nt
        end if

        if (str_endswith(line, 'dihedral types')) then
            read(line,*) nt
            if (.not. allocated(atc%dihedral_specs))  then
                allocate(atc%dihedral_specs(nt))
            else
                allocate(isptmp(dhd_toff+nt))
                isptmp(1:dhd_toff) = atc%dihedral_specs
                call move_alloc(isptmp, atc%dihedral_specs)
            end if

            atc%num_dihedral_types = dhd_toff + nt
        end if

        if (str_endswith(line, 'improper types')) then
            read(line,*) nt
            if (.not. allocated(atc%improper_specs))  then
                allocate(atc%improper_specs(nt))
            else
                allocate(isptmp(imp_toff+nt))
                isptmp(1:imp_toff) = atc%improper_specs
                call move_alloc(isptmp, atc%improper_specs)
            end if

            atc%num_improper_types = imp_toff + nt
        end if

        if (str_endswith(line, 'atoms')) then
            read(line,*) n
            if (.not. allocated(atc%atoms))  then
                allocate(atc%atoms(2,n)); allocate(atc%charge(n))
                allocate(atc%coordinates(3,n))
                allocate(atc%image_box(3,n))
                atc%image_box = 0
            else
                allocate(itmp2(2,atm_idoff+n))
                itmp2(:,1:atm_idoff) = atc%atoms
                call move_alloc(itmp2, atc%atoms)

                allocate(rtmp(atm_idoff+n))
                rtmp(1:atm_idoff) = atc%charge
                call move_alloc(rtmp, atc%charge)

                allocate(rtmp2(3,atm_idoff+n))
                rtmp2(:,1:atm_idoff) = atc%coordinates
                call move_alloc(rtmp2, atc%coordinates)

                allocate(itmp2(3,atm_idoff+n))
                itmp2(:,1:atm_idoff) = atc%image_box
                call move_alloc(itmp2, atc%image_box)
                atc%image_box(:,atm_idoff+1:) = 0
            end if
            atc%num_atoms = atm_idoff + n
        end if

        if (str_endswith(line, 'bonds')) then
            read(line,*) n
            if (.not. allocated(atc%bonds))  then
                allocate(atc%bonds(3,n))
            else
                allocate(itmp2(3,bnd_idoff+n))
                itmp2(:,1:bnd_idoff) = atc%bonds
                call move_alloc(itmp2, atc%bonds)
            end if
            atc%num_bonds = bnd_idoff + n
        end if

        if (str_endswith(line, 'angles')) then
            read(line,*) n
            if (.not. allocated(atc%angles))  then
                allocate(atc%angles(4,n))
            else
                allocate(itmp2(4,ang_idoff+n))
                itmp2(:,1:ang_idoff) = atc%angles
                call move_alloc(itmp2, atc%angles)
            end if
            atc%num_angles = ang_idoff + n
        end if

        if (str_endswith(line, 'dihedrals')) then
            read(line,*) n
            if (.not. allocated(atc%dihedrals))  then
                allocate(atc%dihedrals(5,n))
            else
                allocate(itmp2(5,dhd_idoff+n))
                itmp2(:,1:dhd_idoff) = atc%dihedrals
                call move_alloc(itmp2, atc%dihedrals)
            end if
            atc%num_dihedrals = dhd_idoff + n
        end if

        if (str_endswith(line, 'impropers')) then
            read(line,*) n
            if (.not. allocated(atc%impropers))  then
                allocate(atc%impropers(5,n))
            else
                allocate(itmp2(3,imp_idoff+n))
                itmp2(:,1:imp_idoff) = atc%impropers
                call move_alloc(itmp2, atc%impropers)
            end if
            atc%num_impropers = imp_idoff + n
        end if

        if (str_endswith(line, 'xlo xhi')) read(line,*) xlo, xhi
        if (str_endswith(line, 'ylo yhi')) read(line,*) ylo, yhi
        if (str_endswith(line, 'zlo zhi')) read(line,*) zlo, zhi
        if (str_endswith(line, 'xy xz yz')) then
            read(line,*) xy, xz, yz
            if ((xy /= 0.0_rp) .or. (xz /= 0.0_rp) .or. (yz /= 0.0_rp)) then
                is_triclinic = .true.
            end if
        end if

        !Read body section
        if (line == 'Masses') then
            read(fu, *) !Skip the line following a section header
            do it = atm_toff+1, atc%num_atom_types
                read(fu,*) at_i, atc%atom_specs(at_i+atm_toff)%mass
            end do
        end if

        if (line == 'Atoms') then
            read(fu, *) !Skip the line following a section header
            do it = atm_idoff+1, atc%num_atoms
                read(fu,*) i, imol, typ, atc%charge(atm_idoff+i), &
                    atc%coordinates(:,atm_idoff+i)
                atc%atoms(1,atm_idoff+i) = typ + atm_toff
                if (imol /=0) then
                    atc%atoms(2,atm_idoff+i) = imol + mol_idoff
                else
                    atc%atoms(2,atm_idoff+i) = imol
                end if
                atc%atom_pop(atm_toff+typ) = atc%atom_pop(atm_toff+typ) + 1
            end do
        end if

        if (line == 'Bonds') then
            read(fu, *) !Skip the line following a section header
            do it = bnd_idoff+1, atc%num_bonds
                read(fu,*) i, typ, atc%bonds(2:3,bnd_idoff+i)
                atc%bonds(1,bnd_idoff+i) = typ + bnd_toff
            end do
        end if

        if (line == 'Angles') then
            read(fu, *) !Skip the line following a section header
            do it = ang_idoff+1, atc%num_angles
                read(fu,*) i, typ, atc%angles(2:4,ang_idoff+i)
                atc%angles(1,ang_idoff+i) = typ + ang_toff
            end do
        end if

        if (line == 'Dihedrals') then
            read(fu, *) !Skip the line following a section header
            do it = dhd_idoff+1, atc%num_dihedrals
                read(fu,*) i, typ, atc%dihedrals(2:5,dhd_idoff+i)
                atc%dihedrals(1,dhd_idoff+i) = typ + dhd_toff
            end do
        end if

        if (line == 'Impropers') then
            read(fu, *) !Skip the line following a section header
            do it = imp_idoff+1, atc%num_impropers
                read(fu,*) i, typ, atc%impropers(2:5,imp_idoff+i)
                atc%impropers(1,imp_idoff+i) = typ + imp_toff
            end do
        end if

    end do

    xlo = min(xlo, simbox%xlo); xhi = max(xhi, simbox%xhi)
    ylo = min(ylo, simbox%ylo); yhi = max(yhi, simbox%yhi)
    zlo = min(zlo, simbox%zlo); zhi = max(zhi, simbox%zhi)
    call simbox%init(is_triclinic=is_triclinic)
    call simbox%set_bounds(xlo, xhi, ylo, yhi, zlo, zhi)
    if ( is_triclinic ) call simbox%set_tilt(xy, xz, yz)
    call simbox%update()

    close(fu)

    end subroutine

!******************************************************************************

subroutine write_ldf(simbox, atc, fn_ld, title)
    !! Write to a LAMMPS data file.

    type(smbx_t),     intent(in) :: simbox
    type(atmcfg_t),   intent(in) :: atc
    character(len=*), intent(in) :: fn_ld
        !! Name of the file
    character(len=*), intent(in) :: title
        !! Title of the configuation
    integer :: fu_ld, typ, i

    open(newunit=fu_ld, file=fn_ld, action='write')

    !Header
    write(fu_ld,'(a)') '#'//trim(title)

    write(fu_ld,'(/i0,1x,a)') atc%num_atoms, 'atoms'
    write(fu_ld,'(i0,1x,a)') atc%num_atom_types, 'atom types'

    if (atc%num_bonds > 0) then
        write(fu_ld,'(/i0,1x,a)') atc%num_bonds, 'bonds'
        write(fu_ld,'(i0,1x,a)') atc%num_bond_types, 'bond types'
    end if

    if (atc%num_angles > 0) then
        write(fu_ld,'(/i0,1x,a)') atc%num_angles, 'angles'
        write(fu_ld,'(i0,1x,a)') atc%num_angle_types, 'angle types'
    end if

    if (atc%num_dihedrals > 0) then
        write(fu_ld,'(/i0,1x,a)') atc%num_dihedrals, 'dihedrals'
        write(fu_ld,'(i0,1x,a)') atc%num_dihedral_types, 'dihedral types'
    end if

    if (atc%num_impropers > 0) then
        write(fu_ld,'(/i0,1x,a)') atc%num_impropers, 'impropers'
        write(fu_ld,'(i0,1x,a)') atc%num_improper_types, 'improper types'
    end if

    !Simulation box & tilt factors
    write(fu_ld, *)
    write(fu_ld,'(g0.8,1x,g0.8,1x,a)') simbox%xlo, simbox%xhi, 'xlo xhi'
    write(fu_ld,'(g0.8,1x,g0.8,1x,a)') simbox%ylo, simbox%yhi, 'ylo yhi'
    write(fu_ld,'(g0.8,1x,g0.8,1x,a)') simbox%zlo, simbox%zhi, 'zlo zhi'
    if (simbox%is_triclinic) then
        write(fu_ld,'(3(g0.8,1x),a)') simbox%xy, simbox%xz, simbox%yz, 'xy xz yz'
    end if

    !Body: Masses
    write(fu_ld,*)
    write(fu_ld,'(a)') 'Masses'
    write(fu_ld,*)
    do i = 1, atc%num_atom_types
        write(fu_ld, '(i0,2x,g0.8)') i, atc%atom_specs(i)%mass
    end do

    !Body: Pair Coeffs
!   if (atc%num_vdw_types > 0) then
!       write(fu_ld,*)
!       write(fu_ld,'(a)') 'Pair Coeffs'
!       write(fu_ld,*)
!       !Only interactions between same atom-types
!       do i = 1, atc%num_atom_types
!           typ = i + (2*atc%num_atom_types-i)*(i-1)/2
!           write(fu_ld, '(i0,*(2x,g0.8))') i, atc%vdw_specs(typ)%params
!       end do
!   end if

!   !Body: Bond Coeffs
!   if (atc%num_bond_types > 0) then
!       write(fu_ld,*)
!       write(fu_ld,'(a)') 'Bond Coeffs'
!       write(fu_ld,*)
!       do i = 1, atc%num_bond_types
!           write(fu_ld, '(i0,*(2x,g0.8))') i, atc%bond_specs(i)%params
!       end do
!   end if

!   !Body: Angle Coeffs
!   if (atc%num_angle_types > 0) then
!       write(fu_ld,*)
!       write(fu_ld,'(a)') 'Angle Coeffs'
!       write(fu_ld,*)
!       do i = 1, atc%num_angle_types
!           write(fu_ld, '(i0,*(2x,g0.8))') i, atc%angle_specs(i)%params
!       end do
!   end if

!   !Body: Dihedral Coeffs
!   if (atc%num_dihedral_types > 0) then
!       write(fu_ld,*)
!       write(fu_ld,'(a)') 'Dihedral Coeffs'
!       write(fu_ld,*)
!       do i = 1, atc%num_dihedral_types
!           write(fu_ld, '(i0,*(2x,g0.8))') i, atc%dihedral_specs(i)%params
!       end do
!   end if

!   !Body: Improper Coeffs
!   if (atc%num_improper_types > 0) then
!       write(fu_ld,*)
!       write(fu_ld,'(a)') 'Improper Coeffs'
!       write(fu_ld,*)
!       do i = 1, atc%num_improper_types
!           write(fu_ld, '(i0,*(2x,g0.8))') i, atc%improper_specs(i)%params
!       end do
!   end if

    !Body: Atoms
    write(fu_ld,*)
    write(fu_ld,'(a)') 'Atoms # full'
    write(fu_ld,*)

    if (allocated(atc%image_box)) then
        do i = 1, atc%num_atoms
            write(fu_ld,'(i0,2(1x,i0),4(1x,es13.6),3(1x,i0))')   &
                i, atc%atoms(2,i), atc%atoms(1,i), atc%charge(i), &
                atc%coordinates(:,i), atc%image_box(:,i)
        end do
    else
        do i = 1, atc%num_atoms
            write(fu_ld,'(i0,2(1x,i0),4(1x,es13.6))') i, atc%atoms(2,i), &
                atc%atoms(1,i), atc%charge(i), atc%coordinates(:,i)
        end do
    end if

    !Body: Bonds
    if (atc%num_bonds > 0) then
        write(fu_ld,*)
        write(fu_ld,'(a)') 'Bonds'
        write(fu_ld,*)
        do i = 1, atc%num_bonds
            write(fu_ld,'(i0,2x,3(i0,2x))') i, atc%bonds(:,i)
        end do
    end if

    !Body: Angles
    if (atc%num_angles > 0) then
        write(fu_ld,*)
        write(fu_ld,'(a)') 'Angles'
        write(fu_ld,*)
        do i = 1, atc%num_angles
            write(fu_ld,'(i0,2x,4(i0,2x))') i, atc%angles(:,i)
        end do
    end if

    !Body: Dihedrals
    if (atc%num_dihedrals > 0) then
        write(fu_ld,*)
        write(fu_ld,'(a)') 'Dihedrals'
        write(fu_ld,*)
        do i = 1, atc%num_dihedrals
            write(fu_ld,'(i0,2x,5(i0,2x))') i, atc%dihedrals(:,i)
        end do
    end if

    !Body: Impropers
    if (atc%num_impropers > 0) then
        write(fu_ld,*)
        write(fu_ld,'(a)') 'Impropers'
        write(fu_ld,*)
        do i = 1, atc%num_impropers
            write(fu_ld,'(i0,2x,5(i0,2x))') i, atc%impropers(:,i)
        end do
    end if

    close(fu_ld)

    end subroutine

!*******************************************************************************

!subroutine read_dump(simbox, atc, fn, mxrdln, nts)
!    !! Reads a LAMMPS ascii dump file file. The dump file must contain data for
!    !! for only one timestep.
!
!    type(smbx_t),     intent(out) :: simbox 
!    type(atmcfg_t),   intent(in out) :: atc
!    character(len=*), intent(in)  :: fn
!    integer,          intent(in)  :: mxrdln
!    integer(ip_long), intent(out) :: nts
!    character(len=:), allocatable :: line
!    character(len=:), dimension(:), allocatable :: words, columns
!    integer :: i, ios, icol_atm_id, iatm, ncols, natoms, fu
!    real(rp) :: xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz
!
!    open(newunit=fu, file=fn, action='read', status='old')
!
!    do 
!        call readline(fu, line, '', ios)
!        if (ios /= 0) exit
!
!        if (str_startswith(line, 'ITEM: TIMESTEP')) read(fu,*) nts
!
!        if (str_startswith(line, 'ITEM: NUMBER OF ATOMS')) read(fu,*) natoms
!
!        if (str_startswith(line, 'ITEM: BOX BOUNDS')) then
!            call str_split_to_array(line, words)
!            call simbox%set_boundary_style(trim(words(4)), &
!                trim(words(5)), trim(words(6)))
!            if (size(words)==6) then
!                read(fu,*) xlo, xhi
!                read(fu,*) ylo, yhi
!                read(fu,*) zlo, zhi
!                call simbox%set_bounds(xlo, xhi, ylo, yhi, zlo, zhi)
!            else
!                read(fu,*) xlo, xhi, xy
!                read(fu,*) ylo, yhi, xz
!                read(fu,*) zlo, zhi, yz
!                call simbox%set_bounds(xlo, xhi, ylo, yhi, zlo, zhi)
!                call simbox%set_tilt(xy, xz, yz)
!        end if
!
!        if (str_startswith(line, 'ITEM: ATOMS')) then
!            call str_split_to_array(line, words)
!            allocate(columns, source=words(3:))
!            ncols = size(columns)
!            do icol = 1, ncols
!                if(words(icol)=='id') then
!                    icol_atom_id = icol
!                else if (words(icol)=='x') then
!                    atc%coord_desc(1) = 'uw'
!                else if (words(icol)=='xs') then
!                    atc%coord_desc(1) = 'sw'
!                else if (words(icol)=='xu') then
!                    atc%coord_desc(1) = 'uu'
!                else if (words(icol)=='xsu') then
!                    atc%coord_desc(1) = 'su'
!                else if (words(icol)=='y') then
!                    atc%coord_desc(2) = 'uw'
!                else if (words(icol)=='ys') then
!                    atc%coord_desc(2) = 'sw'
!                else if (words(icol)=='yu') then
!                    atc%coord_desc(2) = 'uu'
!                else if (words(icol)=='ysu') then
!                    atc%coord_desc(2) = 'su'
!                else if (words(icol)=='z') then
!                    atc%coord_desc(3) = 'uw'
!                else if (words(icol)=='zs') then
!                    atc%coord_desc(3) = 'sw'
!                else if (words(icol)=='zu') then
!                    atc%coord_desc(3) = 'uu'
!                else if (words(icol)=='zsu') then
!                    atc%coord_desc(3) = 'su'
!                end if
!            end do
!
!            do i = 1, natoms
!                call readline(fu, line, '', ios)
!                call str_split_to_array(line, words)
!                read(words(icol_atm_id),*) iatm
!                do icol = 1, ncols
!                    select case(columns(icol))
!                    case('q')
!                        read(words(icol),*) atc%charge(iatm)
!                    case('x', 'xs', 'xu', 'xsu')
!                        read(words(icol),*) atc%coordinates(1,iatm)
!                    case('y', 'ys', 'yu', 'ysu')
!                        read(words(icol),*) atc%coordinates(2,iatm)
!                    case('z', 'zs', 'zu', 'zsu')
!                        read(words(icol),*) atc%coordinates(3,iatm)
!                    case('ix')
!                        read(words(icol),*) atc%image_box(1,iatm)
!                    case('iy')
!                        read(words(icol),*) atc%image_box(2,iatm)
!                    case('iz')
!                        read(words(icol),*) atc%image_box(3,iatm)
!                    case('vx')
!                        read(words(icol),*) atc%velocities(1,iatm)
!                    case('vy')
!                        read(words(icol),*) atc%velocities(2,iatm)
!                    case('vz')
!                        read(words(icol),*) atc%velocities(3,iatm)
!                    case('fx')
!                        read(words(icol),*) atc%forces(1,iatm)
!                    case('fy')
!                        read(words(icol),*) atc%forces(2,iatm)
!                    case('fz')
!                        read(words(icol),*) atc%forces(3,iatm)
!                    case('sxx')
!                        read(words(icol),*) atc%stresses(1,iatm)
!                    case('syy')
!                        read(words(icol),*) atc%stresses(2,iatm)
!                    case('szz')
!                        read(words(icol),*) atc%stresses(3,iatm)
!                    case('sxy')
!                        read(words(icol),*) atc%stresses(4,iatm)
!                    case('sxz')
!                        read(words(icol),*) atc%stresses(5,iatm)
!                    case('syz')
!                        read(words(icol),*) atc%stresses(6,iatm)
!                    case default
!                        cycle
!                    end select
!                end do
!            end do
!        end if
!    end do
!
!    close(fu)
!
!    end subroutine

!*******************************************************************************

subroutine write_dmp(simbox, atc, nts, fn)
    !! Writes a single time step to a LAMMPS ascii dump file. 

     type(smbx_t),     intent(in) :: simbox 
     type(atmcfg_t),   intent(in) :: atc
     integer(ip_long), intent(out) :: nts
     character(len=*), intent(in)  :: fn
     integer :: iatm, fu

     open(newunit=fu, file=fn, action='write', status='replace')
     write(fu,*) 'ITEM: TIMESTEP'
     write(fu,*) nts
     write(fu,*) 'ITEM: NUMBER OF ATOMS'
     write(fu,*) atc%num_atoms
     if (simbox%is_triclinic) then
        write(fu,*) 'ITEM: BOX BOUNDS', simbox%boundary
        write(fu,*) simbox%xlo, simbox%xhi
        write(fu,*) simbox%ylo, simbox%yhi
        write(fu,*) simbox%zlo, simbox%zhi
    else
        write(fu,*) 'ITEM: BOX BOUNDS', simbox%boundary, 'xy xz yz'
        write(fu,*) simbox%xlo, simbox%xhi, simbox%xy
        write(fu,*) simbox%ylo, simbox%yhi, simbox%xz
        write(fu,*) simbox%zlo, simbox%zhi, simbox%yz
    end if
    write(fu,*) 'ITEM: ATOMS id type x y z'
    do iatm = 1, atc%num_atoms
        write(fu,*) iatm, atc%atoms(1,iatm), atc%coordinates(:,iatm)
    end do

    close(fu)

    end subroutine

!*******************************************************************************

subroutine write_xyz(atc, fn_xyz, title)
    !! Write to an XYZ file.

    type(atmcfg_t),   intent(in) :: atc
    character(len=*), intent(in) :: fn_xyz
        !! Name of the XYZ file
    character(len=*), intent(in) :: title
        !! Title (for the configuration)
    integer :: iatm, na, at, fu_xyz
    character(len=:), allocatable :: atm_nam

    na = atc%num_atoms

    open(newunit=fu_xyz, file=fn_xyz, action='write', status='replace')

    write(fu_xyz, '(i0)') na
    write(fu_xyz, '(a)') title

    do iatm = 1, na
        at = atc%atoms(1,iatm)
        atm_nam = atc%atom_specs(at)%name
        write(fu_xyz,'(a,3(1x,es22.14))') trim(atm_nam), atc%coordinates(:,iatm)
    end do

    close(fu_xyz)

    end subroutine

!*******************************************************************************

subroutine read_xyz(atc, fn_xyz)
    !! Reads atom configuration from an XYZ file.

    type(atmcfg_t),   intent(in out) :: atc
    character(len=*), intent(in) :: fn_xyz
        !! Name of the XYZ file
    real(rp), dimension(:,:), allocatable :: rtmp2
    integer :: iatm, i, na, fu_xyz, atm_idoff
    character(len=8) :: atm_nam

    open(newunit=fu_xyz, file=fn_xyz, action='read', status='old')

    read(fu_xyz, *) na
    read(fu_xyz, *) !Skip the title line

    atm_idoff = atc%num_atoms

    if (.not. allocated(atc%coordinates))  then
        allocate(atc%coordinates(3,na))
    else
        allocate(rtmp2(3,atm_idoff+na))
        rtmp2(:,1:atm_idoff) = atc%coordinates
        call move_alloc(rtmp2, atc%coordinates)
    end if
    atc%num_atoms = atm_idoff + na

    do i = 1, na
        iatm = i + atm_idoff
        read(fu_xyz,*) atm_nam, atc%coordinates(:,iatm)
    end do

    close(fu_xyz)

    end subroutine

!*******************************************************************************

subroutine read_molrec(fn_mr, molecules)
    !! Reads molecule record
     
    character(len=*), intent(in) :: fn_mr
        !! Name of the molrec file
    integer, dimension(:,:), allocatable, intent(in out) :: molecules
    integer :: fu_mr, i, imol, n, rows, cols

    open(newunit=fu_mr, file=fn_mr, action='read', status='old')
    read(fu_mr,*) n
    if (.not. allocated(molecules)) then
        allocate(molecules(9,n))
    else
        rows = size(molecules,1); cols = size(molecules,2)
        if ((rows/=9) .or. (cols/=n)) then
            deallocate(molecules)
            allocate(molecules(9,n))
        end if
    end if

    !6: bndgrp pop, 7-9: bndgrp
    read(fu_mr,*)
    do i = 1, n
        read(fu_mr, *) imol, molecules(:,i)
    end do

    close(fu_mr)
    
    end subroutine

!*******************************************************************************

subroutine read_molinfo(fn, mol_names, mol_bnds)
    !! Reads molecule records
     
    character(*), intent(in) :: fn
        !! Name of the file containing molecule names and bounds of atom ids
    character(:), dimension(:), allocatable, intent(out) :: mol_names
    integer, dimension(:,:), allocatable, intent(out) :: mol_bnds
    character(:), dimension(:), allocatable :: words
    character(:), allocatable :: line
    integer :: i, ios, imol, fu, mxlen, num_mols

    mxlen = 0
    open(newunit=fu, file=fn, action='read', status='old')
    read(fu,*) !Skip the first line
    do
        call readline(fu, line, '', ios)
        if (ios /= 0) exit

        if (str_startswith(line, 'MOLECULES')) then
            call str_split_to_array(line, words)
            read(words(2),*) num_mols
            allocate(mol_bnds(2,num_mols))
            do i = 1, num_mols
                call readline(fu, line, '', ios)
                call str_split_to_array(line, words)
                mxlen = max( mxlen, len_trim(words(2)) )
            end do
            allocate(character(len=mxlen):: mol_names(num_mols))
        end if
    end do
    close(fu)
    
    open(newunit=fu, file=fn, action='read', status='old')
    read(fu,*) !Skip the first line
    do
        call readline(fu, line, '', ios)
        if (ios /= 0) exit

        if (str_startswith(line, 'MOLECULES')) then
            do i = 1, num_mols
                read(fu,*) imol, mol_names(i), mol_bnds(:,i)
            end do
        end if
    end do
    close(fu)

    end subroutine

!*******************************************************************************

end module config_io_m
