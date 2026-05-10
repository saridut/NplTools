!! Analysis routines for configurations of ligands grafted on nanoplatelets.

module ligang_calculator_m

use constants_m
use utils_math_m
use atmcfg_m
use simbox_m

implicit none

private
public :: ligang_calculator_t

type ligang_calculator_t
    logical  :: is_slab = .false.
    logical  :: normalize = .true.
    real(rp) :: zads = 0.0_rp !Adsorption distance
    real(rp) :: bin_size_theta = 0.0_rp
    real(rp) :: bin_size_phi = 0.0_rp
    real(rp), dimension(2,3) :: xtal_xtnt = 0.0_rp
    real(rp), dimension(:,:), allocatable :: hist_theta
    real(rp), dimension(:,:), allocatable :: hist_phi
    character(:), allocatable :: fn_out
    contains
        procedure :: init => lgac_init
        procedure :: delete => lgac_delete
        procedure :: writeout => lgac_writeout
        procedure :: process_frame => lgac_process_frame
end type ligang_calculator_t


contains

!*******************************************************************************

subroutine lgac_init(this, atc)

    class(ligang_calculator_t), intent(out) :: this
    type(atmcfg_t), intent(in) :: atc
    integer  :: i, at, num_bins, na_xtal
    integer, allocatable :: xtal_atom_types(:), selected_atom_types(:)

    xtal_atom_types = [1,2]
    this%is_slab = .false.
    this%zads = 2.5_rp
    this%bin_size_theta = 3.0_rp
    this%bin_size_phi = 3.0_rp
    this%fn_out = "ang"

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

    !For angle theta
    num_bins = ceiling(math_pi_2/this%bin_size_theta)
    this%bin_size_theta = math_pi_2/num_bins 
    allocate(this%hist_theta(3,num_bins))
    this%hist_theta = 0.0_rp
    do i = 1, num_bins
        this%hist_theta(1,i) = (real(i)-0.5_rp)*this%bin_size_theta
    end do

    !For angle phi
    num_bins = ceiling(math_pi/this%bin_size_phi)
    this%bin_size_phi = math_pi/num_bins 
    allocate(this%hist_phi(3,num_bins))
    this%hist_phi = 0.0_rp
    do i = 1, num_bins
        this%hist_phi(1,i) = (real(i)-0.5_rp)*this%bin_size_phi - math_pi_2 
    end do
    
    end subroutine

!*******************************************************************************

subroutine lgac_delete(this)

    class(ligang_calculator_t), intent(in out) :: this
    
    this%is_slab = .false.
    this%zads = 0.0_rp
    this%bin_size_theta = 0.0_rp
    this%bin_size_phi = 0.0_rp
    this%xtal_xtnt = 0.0_rp

    if (allocated(this%hist_theta)) deallocate(this%hist_theta)
    if (allocated(this%hist_phi)) deallocate(this%hist_phi)

    end subroutine

!*******************************************************************************

subroutine lgac_writeout(this)

    class(ligang_calculator_t), intent(in out) :: this
    integer :: i, fu

    open(newunit=fu, file=trim(this%fn_out)//'_theta.csv', &
        action='write', status='replace')
    write(fu, '(a16,",",1x,a16,",",1x,a16)') 'theta', 'top', 'bot'
    if (this%normalize) then
        this%hist_theta(2,:) = this%hist_theta(2,:) &
                /(this%bin_size_theta*sum(this%hist_theta(2,:)))
        this%hist_theta(3,:) = this%hist_theta(3,:) &
                /(this%bin_size_theta*sum(this%hist_theta(3,:)))
    end if
    do i = 1, size(this%hist_theta,2)
        write(fu, '(3(es16.7,",",1x))') this%hist_theta(:,i)
    end do
    close(fu)

    open(newunit=fu, file=trim(this%fn_out)//'_phi.csv', &
        action='write', status='replace')
    write(fu, '(a16,",",1x,a16,",",1x,a16)') 'phi', 'top', 'bot'
    if (this%normalize) then
        this%hist_phi(2,:) = this%hist_phi(2,:) &
                /(this%bin_size_phi*sum(this%hist_phi(2,:)))
        this%hist_phi(3,:) = this%hist_phi(3,:) &
                /(this%bin_size_phi*sum(this%hist_phi(3,:)))
    end if
    do i = 1, size(this%hist_phi,2)
        write(fu, '(3(es16.7,",",1x))') this%hist_phi(:,i)
    end do
    close(fu)

    end subroutine
    
!*******************************************************************************

subroutine lgac_process_frame(this, simbox, atc, mol_names, mol_bnds)

    class(ligang_calculator_t), intent(in out) :: this
    type(smbx_t), intent(in) :: simbox
    type(atmcfg_t), intent(in) :: atc
    character(*), intent(in) :: mol_names(:)
    integer, intent(in) :: mol_bnds(:,:)
    real(rp), dimension(3) :: rh, rt, dr
    real(rp) :: zi_mag, theta, phi
    real(rp) :: xtal_xlo, xtal_xhi, xtal_ylo, xtal_yhi, xtal_zlo, xtal_zhi
    integer :: imol, iatm, iatm_beg, iatm_head, iatm_tail
    integer :: i, k, at, mt
    logical :: is_top

!   xtal_xlo = this%xtal_xtnt(1,1); xtal_xhi = this%xtal_xtnt(2,1)
!   xtal_ylo = this%xtal_xtnt(1,2); xtal_yhi = this%xtal_xtnt(2,2)
!   xtal_zlo = this%xtal_xtnt(1,3); xtal_zhi = this%xtal_xtnt(2,3)

!   do imol = 1, size(mol_bnds,2)
!       mt = molecules(1, imol)
!       if (mt /= this%mol_type) cycle
!       iatm_beg = molecules(2,imol)
!       iatm_head = iatm_beg + molecules(4, imol) - 1
!       rh = atc%coordinates(:,iatm_head)
!   
!       if (.not. this%is_slab) then
!           if ( (rh(1) < xtal_xlo) .or. (rh(1) > xtal_xhi) ) cycle
!           if ( (rh(2) < xtal_ylo) .or. (rh(2) > xtal_yhi) ) cycle
!       end if

!       if ( rh(3) >= xtal_zhi ) then
!           zi_mag = rh(3) - xtal_zhi
!           is_top = .true.
!       else if ( rh(3) <= xtal_zlo ) then
!           zi_mag = xtal_zlo - rh(3)
!           is_top = .false.
!       else
!           cycle
!       end if
!       if (zi_mag > this%zads) cycle

!       iatm_tail = iatm_beg + molecules(5, imol) - 1
!       rt = atc%coordinates(:,iatm_tail)
!       dr = rt - rh; call simbox%get_image(dr); dr = dr/norm2(dr)
!       theta = acos(dr(3)); phi = atan(dr(2)/dr(1))

!       if (is_top) then
!           if ( (theta >= 0.0_rp) .and. (theta <= math_pi_2) ) then
!               k = floor(theta/this%bin_size_theta) + 1
!               this%hist_theta(2,k) = this%hist_theta(2,k) + 1
!           end if

!           k = floor((phi+math_pi_2)/this%bin_size_phi) + 1
!           this%hist_phi(2,k) = this%hist_phi(2,k) + 1
!       else
!           theta = math_pi - theta
!           if ( (theta >= 0.0_rp) .and. (theta <= math_pi_2) ) then
!               k = floor(theta/this%bin_size_theta) + 1
!               this%hist_theta(3,k) = this%hist_theta(3,k) + 1
!           end if

!           k = floor((phi+math_pi_2)/this%bin_size_phi) + 1
!           this%hist_phi(3,k) = this%hist_phi(3,k) + 1
!       end if
!   end do

    end subroutine

!*******************************************************************************

end module ligang_calculator_m
