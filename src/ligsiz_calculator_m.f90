!! Analysis routines for calculating the distribution of atoms of ligands 
!! grafted on nanoplatelets perpendicular to the NPL.

module zdf_calculator_m

use strings_m
use constants_m
use vector_m
use utils_math_m
use atmcfg_m
use simbox_m

implicit none

type zdf_calculator_t
    logical  :: is_slab = .false.
    logical  :: normalize = .true.
    logical  :: mass_density = .false.
    integer  :: num_frames = 0
    real(rp) :: zmax = 0.0_rp
    real(rp) :: zads = 0.0_rp !Adsorption distance
    real(rp) :: bin_size = 0.0_rp
    real(rp) :: xtal_xtnt(2,3) = 0.0_rp
    type(ivector_t) :: selected_atoms
    real(rp), allocatable :: hist(:,:)
    character(:), allocatable :: fn_out
    contains
        procedure :: init => zdfc_init
        procedure :: delete => zdfc_delete
        procedure :: writeout => zdfc_writeout
        procedure :: process_frame => zdfc_process_frame
end type zdf_calculator_t

contains

!*******************************************************************************

subroutine zdfc_init(this, atc)

    class(zdf_calculator_t), intent(out) :: this
    type(atmcfg_t), intent(in) :: atc
    integer  :: i, at, num_bins, na_xtal
    integer, allocatable :: xtal_atom_types(:), selected_atom_types(:)

    selected_atom_types = [3,4,5,6,7]
    xtal_atom_types = [1,2]
    this%is_slab = .false.
    this%normalize = .false.
    this%mass_density = .true.
    this%num_frames = 0
    this%zmax = 30.0_rp
    this%zads = 2.5_rp
    this%bin_size = 0.05_rp
    this%fn_out = "gz.csv"
        
    call ivector_init(this%selected_atoms)
    do i = 1, atc%num_atoms
        at = atc%atoms(1,i)
        if ( any(selected_atom_types==at) ) call this%selected_atoms%append(i)
    end do
    call this%selected_atoms%shrink_to_fit()

    this%xtal_xtnt(1,:) = huge(1.0_rp); this%xtal_xtnt(2,:) = -huge(1.0_rp)
    na_xtal = 0
    do i = 1, atc%num_atoms
        at = atc%atoms(1,i)
        if ( any(xtal_atom_types==at) ) then
            this%xtal_xtnt(1,1) = min(this%xtal_xtnt(1,1), atc%coordinates(1,i))
            this%xtal_xtnt(1,2) = min(this%xtal_xtnt(1,2), atc%coordinates(2,i))
            this%xtal_xtnt(1,3) = min(this%xtal_xtnt(1,3), atc%coordinates(3,i))

            this%xtal_xtnt(2,1) = max(this%xtal_xtnt(2,1), atc%coordinates(1,i))
            this%xtal_xtnt(2,2) = max(this%xtal_xtnt(2,2), atc%coordinates(2,i))
            this%xtal_xtnt(2,3) = max(this%xtal_xtnt(2,3), atc%coordinates(3,i))
            na_xtal = na_xtal + 1
        end if
    end do
    
    write(*,'(/,a,2x,i0)')   'num_atms_xtal', na_xtal
    write(*,'(a,2x,g0.8,1x,g0.8)') 'xtal_x', this%xtal_xtnt(:,1)
    write(*,'(a,2x,g0.8,1x,g0.8)') 'xtal_y', this%xtal_xtnt(:,2)
    write(*,'(a,2x,g0.8,1x,g0.8)') 'xtal_z', this%xtal_xtnt(:,3)
    
    num_bins = ceiling(this%zmax/this%bin_size)
    !Fix bin_size to make all bins equispaced
    this%bin_size = this%zmax/num_bins 
    
    allocate(this%hist(7,num_bins))
    this%hist = 0.0_rp
    
    do i = 1, num_bins
        this%hist(1,i) = (real(i)-0.5_rp)*this%bin_size
    end do

    write(*,'(a,2x,l2)')   'is_slab', this%is_slab
    write(*,'(a,2x,l2)')   'normalize', this%normalize
    write(*,'(a,2x,l2)')   'mass_density', this%mass_density
    write(*,'(a,2x,g0.8)') 'zmax', this%zmax
    write(*,'(a,2x,g0.8)') 'bin_size', this%bin_size
    write(*,'(a,2x,g0.8)') 'num_bins', num_bins
    write(*,'(a,2x,a,/)')  'fn_out', this%fn_out
    
    end subroutine

!*******************************************************************************

subroutine zdfc_delete(this)

    class(zdf_calculator_t), intent(in out) :: this
    
    this%is_slab = .false.
    this%normalize = .false.
    this%mass_density = .false.
    this%num_frames = 0
    this%zmax = 0.0_rp
    this%zads = 0.0_rp
    this%bin_size = 0.0_rp
    this%xtal_xtnt = 0.0_rp

    call this%selected_atoms%delete()
    if (allocated(this%hist)) deallocate(this%hist)
    if (allocated(this%fn_out)) deallocate(this%fn_out)

    end subroutine

!*******************************************************************************

subroutine zdfc_writeout(this)

    class(zdf_calculator_t), intent(in out) :: this
    integer :: i, fu

    open(newunit=fu, file=this%fn_out, action='write', status='replace')
    write(fu, '(a16,*(",",1x,a16))') 'z', 'gz_top-C1', 'gz_bot-C1',&
        'gz_top-O', 'gz_bot-O', 'gz_top-tail', 'gz_bot-tail'

    this%hist(2,:) = this%hist(2,:)/this%num_frames
    this%hist(3,:) = this%hist(3,:)/this%num_frames
    this%hist(4,:) = this%hist(4,:)/this%num_frames
    this%hist(5,:) = this%hist(5,:)/this%num_frames
    this%hist(6,:) = this%hist(6,:)/this%num_frames
    this%hist(7,:) = this%hist(7,:)/this%num_frames

    if (this%normalize) then
        this%hist(2,:) = this%hist(2,:)/(this%bin_size*sum(this%hist(2,:)))
        this%hist(3,:) = this%hist(3,:)/(this%bin_size*sum(this%hist(3,:)))
    end if

    do i = 1, size(this%hist,2)
        write(fu, '(g0.6, *(",",1x,g0.6))') this%hist(:,i)
    end do

    close(fu)

    end subroutine
    
!*******************************************************************************

subroutine zdfc_process_frame(this, atc, mol_names, mol_bnds)

    class(zdf_calculator_t), intent(in out) :: this
    type(atmcfg_t),          intent(in) :: atc
    character(*),            intent(in) :: mol_names(:)
    integer,                 intent(in) :: mol_bnds(:,:)
    real(rp) :: ri(3), rh(3), zi_mag, rix, riy, riz, rjz, volume, mass, &
                xtal_xlo, xtal_xhi, xtal_ylo, xtal_yhi, xtal_zlo, xtal_zhi
    integer :: imol, iatm, iatm_beg, iatm_head, i, k, at

    xtal_xlo = this%xtal_xtnt(1,1); xtal_xhi = this%xtal_xtnt(2,1)
    xtal_ylo = this%xtal_xtnt(1,2); xtal_yhi = this%xtal_xtnt(2,2)
    xtal_zlo = this%xtal_xtnt(1,3); xtal_zhi = this%xtal_xtnt(2,3)

    volume = (xtal_xhi-xtal_xlo)*(xtal_yhi-xtal_ylo)*this%bin_size

    do i = 1, this%selected_atoms%len
        iatm = this%selected_atoms%get_val(i)
        at = atc%atoms(1,iatm)
        imol = atc%atoms(2,iatm)
        mass = atc%atom_specs(at)%mass
        ri = atc%coordinates(:,iatm)
        rix = ri(1); riy = ri(2); riz = ri(3)

        if (.not. this%is_slab) then
            if ( (rix < xtal_xlo) .or. (rix > xtal_xhi) ) cycle
            if ( (riy < xtal_ylo) .or. (riy > xtal_yhi) ) cycle
        end if
        if ( riz >= xtal_zhi ) then
            zi_mag = riz - xtal_zhi
            if (zi_mag > this%zmax) cycle
            k = floor(zi_mag/this%bin_size) + 1
            if (this%mass_density) then
                if (at==5)  then
                    this%hist(2,k) = this%hist(2,k) + (mass/volume)
                else if (at==6) then
                    this%hist(4,k) = this%hist(4,k) + (mass/volume)
                else
                    this%hist(6,k) = this%hist(6,k) + (mass/volume)
                end if
            else
                this%hist(2,k) = this%hist(2,k) + 1
            end if
        else if ( riz <= xtal_zlo ) then
            zi_mag = xtal_zlo - riz
            if (zi_mag > this%zmax) cycle
            k = floor(zi_mag/this%bin_size) + 1
            if (this%mass_density) then
                if (at==5)  then
                    this%hist(3,k) = this%hist(3,k) + (mass/volume)
                else if (at==6) then
                    this%hist(5,k) = this%hist(5,k) + (mass/volume)
                else
                    this%hist(7,k) = this%hist(7,k) + (mass/volume)
                end if
            else
                this%hist(3,k) = this%hist(3,k) + 1
            end if
        else
            cycle
        end if
    end do

    this%num_frames = this%num_frames + 1

!   do imol = 1, size(mol_bnds,2)
!       iatm_beg = mol_bnds(1,imol)
!       iatm_head = iatm_beg + mol_bnds(4, imol) - 1
!       rh = atc%coordinates(:,iatm_head)

!       if (.not. this%is_slab) then
!           if ( (rh(1) < xtal_xlo) .or. (rh(1) > xtal_xhi) ) cycle
!           if ( (rh(2) < xtal_ylo) .or. (rh(2) > xtal_yhi) ) cycle
!       end if

!       if ( rh(3) >= xtal_zhi ) then
!           zi_mag = rh(3) - xtal_zhi
!       else if ( rh(3) <= xtal_zlo ) then
!           zi_mag = xtal_zlo - rh(3)
!       else
!           cycle
!       end if

!       if (zi_mag > this%zmax) cycle
!       if (zi_mag <= this%zads) cycle

!       do i = 1, size(this%aids)
!           iatm = iatm_beg + this%aids(i) - 1
!           at = atc%atoms(1,iatm)
!           mass = atc%atom_specs(at)%mass
!           ri = atc%coordinates(:,iatm)
!           rix = ri(1); riy = ri(2); riz = ri(3)

!           if (.not. this%is_slab) then
!               if ( (rix < xtal_xlo) .or. (rix > xtal_xhi) ) cycle
!               if ( (riy < xtal_ylo) .or. (riy > xtal_yhi) ) cycle
!           end if
!           if ( riz >= xtal_zhi ) then
!               zi_mag = riz - xtal_zhi
!               if (zi_mag > this%zmax) cycle
!               k = floor(zi_mag/this%bin_size) + 1
!               if (this%mass_density) then
!                   this%hist(2,k) = this%hist(2,k) + (mass/volume)
!               else
!                   this%hist(2,k) = this%hist(2,k) + 1
!               end if
!           else if ( riz <= xtal_zlo ) then
!               zi_mag = xtal_zlo - riz
!               if (zi_mag > this%zmax) cycle
!               k = floor(zi_mag/this%bin_size) + 1
!               if (this%mass_density) then
!                   this%hist(2,k) = this%hist(2,k) + (mass/volume)
!               else
!                   this%hist(2,k) = this%hist(2,k) + 1
!               end if
!           else
!               cycle
!           end if
!       end do
!   end do

    end subroutine

!*******************************************************************************

end module zdf_calculator_m
