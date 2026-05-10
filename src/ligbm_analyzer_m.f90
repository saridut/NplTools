!! Analysis routines for configurations of ligands grafted on nanoplatelets.

module ligbm_analyzer_m

use constants_m
use utils_math_m
use atmcfg_m
use simbox_m

implicit none

type ligbm_analyzer_t
    logical  :: is_slab = .false.
    logical  :: normalize = .true.
    real(rp) :: zads = 0.0_rp
    real(rp) :: bin_size_theta = 0.0_rp
    real(rp) :: bin_size_phi = 0.0_rp
    real(rp) :: bin_size_psi = 0.0_rp
    real(rp), dimension(2,3) :: xtal_xtnt = 0.0_rp
    real(rp), dimension(:,:), allocatable :: hist_theta
    real(rp), dimension(:,:), allocatable :: hist_phi
    real(rp), dimension(:,:), allocatable :: hist_psi
    character(:), allocatable :: fn_out
    contains
        procedure :: init => lgbma_init
        procedure :: delete => lgbma_delete
        procedure :: writeout => lgbma_writeout
        procedure :: process_frame => lgbma_process_frame
end type ligbm_analyzer_t


contains

!*******************************************************************************

subroutine lgbma_init(this, atc)

    class(ligbm_analyzer_t), intent(out) :: this
    type(atmcfg_t), intent(in) :: atc
    integer  :: i, at, num_bins, na_xtal
    integer, allocatable :: xtal_atom_types(:), selected_atom_types(:)

    this%is_slab = .false.
    this%zads = 2.5_rp
    this%bin_size_theta = 0.05 !in radians
    this%bin_size_phi = 0.05 !in radians
    this%bin_size_psi = 0.05 !in radians
    this%fn_out = "bm"

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

    !For angle psi
    num_bins = ceiling(math_pi/this%bin_size_psi)
    this%bin_size_psi = math_pi/num_bins 
    allocate(this%hist_psi(3,num_bins))
    this%hist_psi = 0.0_rp
    do i = 1, num_bins
        this%hist_psi(1,i) = (real(i)-0.5_rp)*this%bin_size_psi - math_pi_2 
    end do
    
    end subroutine

!*******************************************************************************

subroutine lgbma_delete(this)

    class(ligbm_analyzer_t), intent(in out) :: this
    
    this%is_slab = .false.
    this%zads = 0.0_rp
    this%bin_size_theta = 0.0_rp
    this%bin_size_phi = 0.0_rp
    this%bin_size_psi = 0.0_rp
    this%xtal_xtnt = 0.0_rp

    if (allocated(this%hist_theta)) deallocate(this%hist_theta)
    if (allocated(this%hist_phi)) deallocate(this%hist_phi)
    if (allocated(this%hist_psi)) deallocate(this%hist_psi)

    end subroutine

!*******************************************************************************

subroutine lgbma_writeout(this)

    class(ligbm_analyzer_t), intent(in out) :: this
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

    open(newunit=fu, file=trim(this%fn_out)//'_psi.csv', &
        action='write', status='replace')
    write(fu, '(a16,",",1x,a16,",",1x,a16)') 'psi', 'top', 'bot'
    if (this%normalize) then
        this%hist_psi(2,:) = this%hist_psi(2,:) &
                /(this%bin_size_psi*sum(this%hist_psi(2,:)))
        this%hist_psi(3,:) = this%hist_psi(3,:) &
                /(this%bin_size_psi*sum(this%hist_psi(3,:)))
    end if
    do i = 1, size(this%hist_psi,2)
        write(fu, '(3(es16.7,",",1x))') this%hist_psi(:,i)
    end do
    close(fu)

    end subroutine
    
!*******************************************************************************

subroutine lgbma_process_frame(this, simbox, atc, mol_names, mol_bnds)

    class(ligbm_analyzer_t), intent(in out) :: this
    type(smbx_t),               intent(in) :: simbox
    type(atmcfg_t),             intent(in) :: atc
    character(*), intent(in) :: mol_names(:)
    integer, intent(in) :: mol_bnds(:,:)
    real(rp), dimension(3) :: rh, ri, rj, rk, dr_ij, dr_ik, dr_jk, uhat
    real(rp) :: zi_mag, theta, phi, psi
    real(rp) :: xtal_xlo, xtal_xhi, xtal_ylo, xtal_yhi, xtal_zlo, xtal_zhi
    integer :: imol, iatm, iatm_beg, iatm_head, iatm_tail, iatm_i, iatm_j, iatm_k
    integer :: i, k, at, mt
    logical :: is_top

    xtal_xlo = this%xtal_xtnt(1,1); xtal_xhi = this%xtal_xtnt(2,1)
    xtal_ylo = this%xtal_xtnt(1,2); xtal_yhi = this%xtal_xtnt(2,2)
    xtal_zlo = this%xtal_xtnt(1,3); xtal_zhi = this%xtal_xtnt(2,3)

!   do imol = 1, size(mol_bnds,2)
!       mt = molecules(1, imol)
!       if (mt /= this%mol_type) cycle
!       iatm_beg = mol_bnds(1,imol)
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

!       iatm_i = iatm_beg + molecules(8, imol) - 1 !Central atom
!       iatm_j = iatm_beg + molecules(7, imol) - 1
!       iatm_k = iatm_beg + molecules(9, imol) - 1

!       ri = atc%coordinates(:,iatm_i)
!       rj = atc%coordinates(:,iatm_j)
!       rk = atc%coordinates(:,iatm_k)

!       dr_ij = rj - ri; call simbox%get_image(dr_ij); dr_ij = dr_ij/norm2(dr_ij)
!       dr_ik = rk - ri; call simbox%get_image(dr_ik); dr_ik = dr_ik/norm2(dr_ik)
!       dr_jk = rk - rj; call simbox%get_image(dr_jk); dr_jk = dr_jk/norm2(dr_jk)
!       psi = acos(dot_product(dr_jk, dr_ij)) - math_pi_2
!       call cross(dr_ij, dr_ik, uhat)
!       uhat = uhat/norm2(uhat)
!       theta = acos(uhat(3)); phi = atan(uhat(2)/uhat(1))

!       if (is_top) then
!           if ( (theta >= 0.0_rp) .and. (theta <= math_pi_2) ) then
!               k = floor(theta/this%bin_size_theta) + 1
!               this%hist_theta(2,k) = this%hist_theta(2,k) + 1
!           end if

!           k = floor((phi+math_pi_2)/this%bin_size_phi) + 1
!           this%hist_phi(2,k) = this%hist_phi(2,k) + 1

!           k = floor((psi+math_pi_2)/this%bin_size_psi) + 1
!           this%hist_psi(2,k) = this%hist_psi(2,k) + 1
!       else
!           theta = math_pi - theta
!           if ( (theta >= 0.0_rp) .and. (theta <= math_pi_2) ) then
!               k = floor(theta/this%bin_size_theta) + 1
!               this%hist_theta(3,k) = this%hist_theta(3,k) + 1
!           end if

!           k = floor((phi+math_pi_2)/this%bin_size_phi) + 1
!           this%hist_phi(3,k) = this%hist_phi(3,k) + 1

!           k = floor((psi+math_pi_2)/this%bin_size_psi) + 1
!           this%hist_psi(3,k) = this%hist_psi(3,k) + 1
!       end if
!   end do

    end subroutine

!*******************************************************************************

end module ligbm_analyzer_m
