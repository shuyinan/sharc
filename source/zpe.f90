!******************************************
!
!    SHARC Program Suite
!
!    Copyright (c) 2023 University of Vienna
!
!    This file is part of SHARC.
!
!    SHARC is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published
!    by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHARC is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    inside the SHARC manual.  If not, see
!    <http://www.gnu.org/licenses/>.
!
!******************************************

!> # Module  zpe
!> perform zero point energy correction on semiclassical dynamics
!> this module includes the following algorithms:
!> ZPE pumping - not working yet.
!> LP-ZPE and improved LP-ZPE algorithms
!> 
!> \author Yinan Shu
!> \date Nov. 7, 2023
!>
!> ZPE pumping is implemented on Nov. 18, 2021
!> LP-ZPE  is implemented on Sep. 12, 2022
!> improved LP-ZPE is implemented on Nov.7, 2023

module zpe
  implicit none

  ! ZPE pumping
  public initial_zpe
  public ZPEcorrection
  public ZPEpumping

  private compute_local_mode
  private pumping_state_switch
  private compute_pumping_direction_time
  private propagate_pumping_coeff

  ! LP-ZPE
  public LP_ZPE
  
  private LP_ZPE_MB
  private LP_ZPE_ST
  
  private compute_parallel_velocity
  private correct_bond_velocity
  private optimize_veloc_correction
  private objective_function
  private compute_dlp_zpe

contains

! ===========================================================
! the wrap up for ZPE correction
  subroutine initial_zpe(traj,ctrl)
    use definitions
    use matrix
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl

    integer :: istate

    if (printlevel>3) then
      write(u_log,*) '============================================================='
      write(u_log,*) '    Initialization Zero Point Energy Correction Scheme'
      write(u_log,*) '============================================================='
    endif

    if (ctrl%zpe_correction==1) then ! Do ZPE Pumping
      ! initialize all pumping state to original state, which is state 1
      traj%state_pumping_s=1
 
      ! initialize all coefficients of state 1 to original state, state 2 to zero 
      do istate=1, ctrl%nstates
        traj%coeff_zpe_s2(istate,1)=dcmplx(1.0d0,0.d0)
        traj%coeff_zpe_s2(istate,2)=dcmplx(0.d0,0.d0)
      enddo 
    else if (ctrl%zpe_correction==2) then ! Do LP-ZPE
      traj%lpzpe_cycle(:)=0
      traj%lpzpe_ke_ah(:)=0.d0
      traj%lpzpe_ke_bc(:)=0.d0  
      traj%lpzpe_iter_incycle(:)=1
      traj%lpzpe_iter_outcycle(:)=1
      traj%in_cycle(:)=0
      traj%lpzpe_starttime=0.d0
      traj%lpzpe_endtime=0.d0
    endif

  endsubroutine


! ===========================================================
! the wrap up for ZPE correction
  subroutine ZPEcorrection(traj,ctrl)
    use definitions
    use matrix
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl


    if (ctrl%zpe_correction==1) then ! Do ZPE Pumping
      call ZPEpumping(traj,ctrl)
    else if (ctrl%zpe_correction==2) then ! Do LP-ZPE 
      call LP_ZPE(traj,ctrl)
    endif 

  endsubroutine

! ===========================================================
! Start ZPE pumping algorithm 
! ===========================================================

! ===========================================================
! Using ZPE pumping method to correct each electronic state. 
  subroutine ZPEpumping(traj,ctrl)
    use definitions
    use matrix
    use nuclear 
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl

    complex*16 :: Hpump(ctrl%nstates,2,2), Gradpump(ctrl%nstates,2,ctrl%natom,3)
    real*8 :: svec_pumping_sad(ctrl%nstates,ctrl%natom,3), pumpingtime_s(ctrl%nstates)
    complex*16 :: Rpumping(ctrl%nstates,2,2)
    complex*16 :: denpump(ctrl%nstates,2,2)

    real*8 :: veloc_app_ad(ctrl%natom,3)
    real*8 :: iso_veloc_app_ad(ctrl%natom,3)
    real*8 :: iso_veloc_ad(ctrl%natom,3), iso_geom_ad(ctrl%natom,3)
    real*8 :: iso_veloc_old_ad(ctrl%natom,3)
    real*8 :: iso_accel_ad(ctrl%natom,3)

    integer :: istate, jstate, ipstate, jpstate
    integer :: iatom, idir
    character(8000) :: string

    if (printlevel>3) then
      write(u_log,*) '============================================================='
      write(u_log,*) '                Zero Point Energy Pumping'
      write(u_log,*) '============================================================='
    endif
   
    ! initialize 
    traj%Evib_local_s=dcmplx(0.d0,0.d0)
    traj%Ezpe_local_s=dcmplx(0.d0,0.d0)
    traj%grad_Ezpe_local_sad=0.d0
    
    ! compute forward velocity 
    if (printlevel>4) write(u_log,*) 'Approximate velocity with forward velocity verlet'
    call VelocityVerlet_vstep_approximate(traj%veloc_old_ad, veloc_app_ad, traj%accel_ad, ctrl%dtstep, ctrl%natom)
 
    ! Compute isoinertial coordiantes, velocities and accelerations
    ! Notice we use scaled mass mu=1.0
    do iatom=1,ctrl%natom
      do idir=1,3 
        iso_geom_ad(iatom,idir)=traj%geom_ad(iatom,idir)*sqrt(traj%mass_a(iatom))
        iso_veloc_ad(iatom,idir)=traj%veloc_ad(iatom,idir)*sqrt(traj%mass_a(iatom))
        iso_veloc_old_ad(iatom,idir)=traj%veloc_old_ad(iatom,idir)*sqrt(traj%mass_a(iatom))
        iso_veloc_app_ad(iatom,idir)=veloc_app_ad(iatom,idir)*sqrt(traj%mass_a(iatom))
        iso_accel_ad(iatom,idir)=traj%accel_ad(iatom,idir)*sqrt(traj%mass_a(iatom))
      enddo
    enddo 

    ! Compute local mode energies and gradients. 
    call compute_local_mode(ctrl%natom, ctrl%nstates, iso_veloc_ad, &
      &iso_veloc_old_ad, iso_veloc_app_ad, iso_accel_ad, traj%mass_a, traj%Evib_local_s, &
      &traj%Ezpe_local_s, traj%grad_Ezpe_local_sad, traj%grad_MCH_sad, &
      &traj%grad_MCH_old_sad, traj%grad_MCH_old2_sad, traj%grad_ad, &
      &ctrl%dtstep, ctrl%dtstep_old, traj%step)

    if (printlevel>4) then
      call vecwrite(ctrl%nstates, traj%Evib_local_s, u_log, 'vibrational energy of local mode', 'F14.9')
      call vecwrite(ctrl%nstates, traj%Ezpe_local_s, u_log, 'zero point energy of local mode', 'F14.9')
      do istate=1,ctrl%nstates
        write(string,'(A20,I3)') 'Ezpe Gradient state ',istate
        call vec3write(ctrl%natom, traj%grad_Ezpe_local_sad(istate,:,:), u_log, trim(string), 'F14.9')
      enddo
    endif

    ! Construction pumping Hamiltonian and gradient
    Hpump=dcmplx(0.d0,0.d0)
    Gradpump=dcmplx(0.d0,0.d0)
    do istate=1,ctrl%nstates
      Hpump(istate,1,1)=traj%H_MCH_ss(istate,istate)
      Hpump(istate,2,2)=traj%H_MCH_ss(istate,istate)-traj%Ezpe_local_s(istate)
      Gradpump(istate,1,:,:)=traj%grad_MCH_sad(istate,:,:)
      Gradpump(istate,2,:,:)=traj%grad_MCH_sad(istate,:,:)-traj%grad_Ezpe_local_sad(istate,:,:)
    enddo

    ! Pumping state switch
    call pumping_state_switch(ctrl%nstates, traj%Evib_local_s, traj%Ezpe_local_s, &
      &traj%state_pumping_s)

    if (printlevel>4) then
      write(u_log,*) 'pumping state:'
      write(u_log,*) (traj%state_pumping_s(istate),istate=1,ctrl%nstates)
    endif

    ! Compute pumping direction and time
    call compute_pumping_direction_time(ctrl%natom, ctrl%nstates, traj%veloc_ad,&
      &traj%mass_a, traj%Evib_local_s, traj%Ezpe_local_s, traj%state_pumping_s, &
      &svec_pumping_sad, pumpingtime_s)

    if (printlevel>4) then
      call vecwrite(ctrl%nstates, pumpingtime_s, u_log, 'pumping time', 'F14.9')
    endif

    
    ! propagate pumping population, decay of mixing
    if (printlevel>4) then
      write(u_log,*) 'old pumping coefficients:'
      do istate=1,ctrl%nstates
        write(u_log,'(2(F14.9,1X),2X,2(F14.9,1X))') traj%coeff_zpe_s2(istate,1), traj%coeff_zpe_s2(istate,2)
      enddo
      do istate=1, ctrl%nstates
        do ipstate=1, 2
          do jpstate=1, 2
            denpump(istate,ipstate,jpstate)=traj%coeff_zpe_s2(istate,ipstate)*conjg(traj%coeff_zpe_s2(istate,jpstate))
          enddo
        enddo
      enddo
      write(u_log,*) 'old diagonal pumping density:'
      do istate=1,ctrl%nstates
        write(u_log,'(2(F14.9,1X),2X,2(F14.9,1X))') denpump(istate,1,1), denpump(istate,2,2)
      enddo
    endif

    call propagate_pumping_coeff(ctrl%nstates, ctrl%dtstep, ctrl%nsubsteps, &
      &pumpingtime_s, traj%state_pumping_s, traj%coeff_zpe_s2, Rpumping) 

    ! compute pumping density matrix
    do istate=1, ctrl%nstates
      do ipstate=1, 2
        do jpstate=1, 2
          denpump(istate,ipstate,jpstate)=traj%coeff_zpe_s2(istate,ipstate)*conjg(traj%coeff_zpe_s2(istate,jpstate))
        enddo
      enddo
    enddo

    if (printlevel>4) then
      write(u_log,*) 'new pumping coefficients:'
      do istate=1,ctrl%nstates
        write(u_log,'(2(F14.9,1X),2X,2(F14.9,1X))') traj%coeff_zpe_s2(istate,1), traj%coeff_zpe_s2(istate,2)
      enddo
      do istate=1,ctrl%nstates
        write(u_log,*) 'state: ', istate
        call matwrite(2, Rpumping(istate,:,:), u_log, 'Pumping propagator','F14.9')
      enddo 
      write(u_log,*) 'new diagonal pumping density:'
      do istate=1,ctrl%nstates
        write(u_log,'(2(F14.9,1X),2X,2(F14.9,1X))') denpump(istate,1,1), denpump(istate,2,2)
      enddo
    endif

    ! adjust the true Hamiltonian and gradients by pumping Hamiltonian and gradients
    do istate=1,ctrl%nstates
      traj%H_MCH_ss(istate,istate)=(denpump(istate,1,1)+denpump(istate,2,2))*traj%H_MCH_ss(istate,istate)-&
        &denpump(istate,2,2)*traj%Ezpe_local_s(istate)
      traj%grad_MCH_sad(istate,:,:)=(denpump(istate,1,1)+denpump(istate,2,2))*traj%grad_MCH_sad(istate,:,:)-&
        &denpump(istate,2,2)*traj%grad_Ezpe_local_sad(istate,:,:)
    enddo


  endsubroutine


! ==========================================================
!> this subroutine computes the local mode vibration energy and zpe 
!> compute_local_mode(ctrl%natom, ctrl%nstates, iso_veloc_ad, &
!>      &iso_veloc_old_ad, iso_veloc_app_ad, iso_accel_ad, traj%mass_a, traj%Evib_local_s, &
!>      &traj%Ezpe_local_s, traj%grad_Ezpe_local_sad, traj%grad_MCH_sad, &
!>      &traj%grad_MCH_old_sad, traj%grad_MCH_old2_sad, traj%grad_ad, &
!>      &ctrl%dtstep, ctrl%dtstep_old, traj%step)
  subroutine compute_local_mode(n, ns, v, vold, v_app, a, m, Evib, Ezpe, gEzpe, g, gold, gold2, traj_g, dt, dtold, step)
    use definitions, only: u_log, printlevel
    use matrix
    use nuclear, only: VelocityVerlet_vstep_approximate
    implicit none
    integer, intent(in) :: n, ns
    real*8, intent(in) :: v(n,3), vold(n,3), v_app(n,3), a(n,3), m(n)
    complex*16, intent(out) :: Evib(ns), Ezpe(ns)
    real*8, intent(out) :: gEzpe(ns,n,3)
    real*8, intent(in) :: g(ns,n,3), gold(ns,n,3), gold2(ns,n,3), traj_g(n,3)
    real*8, intent(in) :: dt, dtold
    integer, intent(in) :: step

    real*8 :: gv_old(ns), gv(ns) ! dV/dt
    real*8 :: eh(ns) ! d2V/dt2
    real*8 :: Ekin_app, mass_app
    real*8 :: normv, normg
    real*8 :: vg_app_ad(n,3)
    integer :: istate, iatom, idir



    if (printlevel>3) then
      write(u_log,*) '============================================================='
      write(u_log,*) 'Calculating local vibration mode energy and elevated surface'
      write(u_log,*) '============================================================='
    endif

    ! compute local mode force constant
    if (step>=1) then 
      gv_old=0.d0
      gv=0.d0
      normv=0.d0
      normg=0.d0
 
      ! here the gradient need to be re-scaled to iso coordiantes 
      do iatom=1,n
        do idir=1,3
          gv_old=gv_old+gold(:,iatom,idir)/sqrt(m(iatom))*vold(iatom,idir)
          gv=gv+g(:,iatom,idir)/sqrt(m(iatom))*v_app(iatom,idir)
          normv=normv+v_app(iatom,idir)**2
          normg=normg+(traj_g(iatom,idir)/sqrt(m(iatom)))**2
        enddo
      enddo
      normv=sqrt(normv)
      normg=sqrt(normg)
      eh=2.0d0*(gv-gv_old)/(dt+dtold)
 
      ! compute velocity along gradient direction 
      do iatom=1,n
        do idir=1,3
          vg_app_ad(iatom,idir)=v_app(iatom,idir)*(traj_g(iatom,idir)/sqrt(m(iatom)))/normg
        enddo
      enddo

      if (printlevel>4) then
        call vec3write(n, v_app, u_log, 'original velocity', 'F14.9')
        call vec3write(n, vg_app_ad, u_log, 'velocity along g direction', 'F14.9')
      endif
          
      Ekin_app=0.d0
      do iatom=1,n
        Ekin_app=Ekin_app+0.5d0*1.d0*sum(v_app(iatom,:)**2)
      enddo
      write(u_log,*) 'original kinetic energy', Ekin_app

      ! compute local vibration mode energy
      ! notice here the scaled mass is 1
      Ekin_app=0.d0
      do iatom=1,n
        Ekin_app=Ekin_app+0.5d0*1.d0*sum(vg_app_ad(iatom,:)**2)
      enddo
      write(u_log,*) 'kinetic energy', Ekin_app

      do istate=1,ns
        if (eh(istate)>0.d0) then
          Evib(istate)=Ekin_app+0.5d0*(gv(istate)**2)/eh(istate)
          Ezpe(istate)=0.5d0*sqrt(eh(istate)/1.0)
          ! compute the derivative of local vibration mode energy 
          if (step==1) then 
            gEzpe(istate,:,:)=(g(istate,:,:)-gold(istate,:,:))/dt
          else if (step>=2) then
            gEzpe(istate,:,:)=2.0d0*((g(istate,:,:)-gold(istate,:,:))/dt-(gold(istate,:,:)-gold2(istate,:,:))/dtold)/(dt+dtold)
          endif 
        else 
          Evib(istate)=Ekin_app
          Ezpe(istate)=0.d0
          gEzpe(istate,:,:)=0.d0
      endif
      enddo ! do istate=1,ns

     endif ! if (step>=1) then

  endsubroutine

! ==========================================================
!> this subroutine switches pumping state
  subroutine pumping_state_switch(ns, Evib, Ezpe, k)
    use definitions, only: u_log, printlevel
    implicit none
    integer, intent(in) :: ns
    complex*16, intent(in) :: Evib(ns), Ezpe(ns)
    integer, intent(inout) :: k(ns) 
   
    integer :: istate

    ! the pumping state switches according to relative value of potential and Ezpe
    do istate=1, ns
      if (real(Evib(istate)) .lt. real(Ezpe(istate))) then 
        k(istate)=2
      else if (real(Evib(istate)) .ge. real(Ezpe(istate))) then 
        k(istate)=1
      endif
    enddo 

  endsubroutine

! ==========================================================
!> this subroutine computes pumping direction and time
!> compute_pumping_direction_time(ctrl%natom, ctrl%nstates, traj%veloc_ad,&
!>   &traj%mass_a, traj%Evib_local_s, traj%Ezpe_local_s, traj%state_pumping_s, &
!>   &svec_pumping_sad, pumpingtime_s)
  subroutine compute_pumping_direction_time(n, ns, v, m, Evib, Ezpe, k, s, t) 
    use definitions, only: u_log, printlevel
    implicit none
    integer, intent(in) :: n, ns
    real*8, intent(in) :: v(n,3), m(n)
    complex*16, intent(in) :: Evib(ns), Ezpe(ns)
    integer, intent(in) :: k(ns)
    real*8, intent(out) :: s(ns,n,3)
    real*8, intent(inout) :: t(ns)

    real*8 :: normv
    real*8 :: p(n,3)
    real*8 :: smag, us(ns), es(ns)
    integer :: istate, iatom, idir

    ! compute norm of v 
    normv=0.d0
    do iatom=1,n
      do idir=1,3
        normv=normv+v(iatom,idir)**2
      enddo
    enddo
    normv=sqrt(normv)
 
    smag=0.d0
    do istate=1,ns

      s(istate,:,:)=v/normv
      do iatom=1,n
        do idir=1,3
          smag=smag+(s(istate,iatom,idir)**2)/m(iatom)
        enddo
      enddo

    us(istate)=0.d0
    do iatom=1,n
      do idir=1,3
        us(istate)=us(istate)+s(istate,iatom,idir)*v(iatom,idir)
      enddo
    enddo
    es(istate)=(us(istate)**2)/(2.d0*smag**2)
    t(istate)=dabs(real(Ezpe(istate)))/(1.d0+1.d0/es(istate))

   enddo !do istate=1,ns

  endsubroutine


! ==========================================================
!> this subroutine propagates the pumping coefficients
  subroutine propagate_pumping_coeff(ns, dt, nsub, tau, k, C, R)
    use definitions, only: u_log, printlevel
    use matrix
    implicit none
    integer, intent(in) :: ns
    integer, intent(in) :: nsub
    real*8, intent(in) :: dt
    real*8, intent(in) :: tau(ns)
    integer, intent(in) :: k(ns)
    complex*16, intent(inout) :: C(ns,2)
    complex*16, intent(inout) :: R(ns,2,2)

    integer :: istep
    real*8 :: dtsubstep
    complex*16 :: Rtotal(ns,2,2), Rprod(ns,2,2)
    complex*16 :: DP(ns,2,2), Rexpd(ns,2,2)
    complex*16 :: Ctmp1(ns,2), Cold(ns,2)
    complex*16 :: ii=dcmplx(0.d0,1.d0)
    complex*16 :: Rtmp(2,2)
 
    integer :: istate, ipstate
    integer :: tracker

    !initialize Ctmp1
    Ctmp1=C
    Cold=C

    dtsubstep=dt/nsub
    Rtotal=dcmplx(0.d0,0.d0)
    do istate=1,ns
      do ipstate=1,2
        Rtotal(istate,ipstate,ipstate)=dcmplx(1.d0,0.d0)
      enddo
    enddo 
    
    ! tracker is used to avoid den(k,k) being exact zero 
    ! if this happens, we use a single state decay, and the pointer state
    ! density is computed by 1-den(i,i)
    tracker=0
    do istate=1, ns
      do ipstate=1, 2
        if (Ctmp1(istate,k(istate)) .eq. dcmplx(0.d0,0.d0)) tracker=1
      enddo
    enddo

    if (tracker==0) then ! do two state case  
      do istep=1,nsub
        DP=dcmplx(0.d0,0.d0)
        do istate=1, ns
          do ipstate=1, 2
            if (ipstate .ne. k(istate)) then
              DP(istate,k(istate),k(istate))=DP(istate,k(istate),k(istate))+&
                &(1.d0,0.d0)*tau(istate)*real(Ctmp1(istate,ipstate)*conjg(Ctmp1(istate,ipstate)))
              DP(istate,ipstate,ipstate)=-0.5d0*(1.d0,0.d0)*tau(istate)
            endif
          enddo
          DP(istate,k(istate),k(istate))=0.5d0*DP(istate,k(istate),k(istate))/&
            &(real(Ctmp1(istate,k(istate))*conjg(Ctmp1(istate,k(istate)))))
        enddo    
        Rexpd=dtsubstep*(DP)
        do istate=1, ns    
          do ipstate=1,2
            Rexpd(istate,ipstate,ipstate)=exp(Rexpd(istate,ipstate,ipstate))
          enddo
          !call exponentiate(2, Rexpd(istate,:,:), (1.d0,0.d0))
          call matmultiply(2, Rexpd(istate,:,:), Rtotal(istate,:,:), Rprod(istate,:,:), 'nn')
          Rtotal(istate,:,:)=Rprod(istate,:,:)
          call matvecmultiply(2, Rtotal(istate,:,:), Cold(istate,:), Ctmp1(istate,:), 'n')
        enddo
      enddo ! do istep=1,nsub
      do istate=1, ns
        call matvecmultiply(2, Rtotal(istate,:,:), Cold(istate,:), C(istate,:), 'n')
      enddo 
      R=Rtotal
    else ! do one state decay 
      do istep=1,nsub
        DP=dcmplx(0.d0,0.d0)
        do istate=1, ns
          do ipstate=1, 2
            if (ipstate .ne. k(istate)) then
              DP(istate,ipstate,ipstate)=-0.5d0*(1.d0,0.d0)*tau(istate)
            endif
          enddo
        enddo
        Rexpd=dtsubstep*(DP)
        do istate=1, ns
          do ipstate=1,2
            Rexpd(istate,ipstate,ipstate)=exp(Rexpd(istate,ipstate,ipstate))
          enddo
          call matmultiply(2, Rexpd(istate,:,:), Rtotal(istate,:,:), Rprod(istate,:,:), 'nn')
          Rtotal(istate,:,:)=Rprod(istate,:,:)
          call matvecmultiply(2, Rtotal(istate,:,:), Cold(istate,:), Ctmp1(istate,:), 'n')
        enddo
      enddo
    do istate=1, ns
      call matvecmultiply(2, Rtotal(istate,:,:), Cold(istate,:), C(istate,:), 'n')
      if (k(istate)==1) then
        C(istate,1)=sqrt(1.d0-C(istate,2)*conjg(C(istate,2)))
      else if (k(istate)==2) then
        C(istate,2)=sqrt(1.d0-C(istate,1)*conjg(C(istate,1)))
      endif 
    enddo
    R=Rtotal
   endif ! if (tracker==0) then

  endsubroutine

! ===========================================================
! End ZPE pumping algorithm 
! ===========================================================


! ===========================================================
! Start LP-ZPE algorithms
! ===========================================================
! ==========================================================
!> main driver of LP-ZPE correction 
  subroutine LP_ZPE(traj,ctrl)
    use definitions
    use matrix
    use nuclear, only: Calculate_ekin
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl

    real*8 :: tstart(ctrl%lpzpe_nah), tend(ctrl%lpzpe_nah)

    integer :: atom1, atom2
    real*8 :: parallel_vec(3), parallel_veloc1, parallel_veloc2
    real*8 :: ke_bond_ah
    integer :: ibond, iatom, jatom, katom, idir

    real*8 :: parallel_vec_tmp(3)
    real*8 :: parallel_veloc1_tmp, parallel_veloc2_tmp
    ! bond energy
    real*8 :: bond_Ekin_AH, bond_Ekin_AH_new
    ! kinetic energy 
    real*8 :: Ekin, Ekin_new
    ! momentum
    real*8 :: momentum(ctrl%natom,3), momentum_new(ctrl%natom,3)
    ! linear momentum, i.e. center of mass momentum 
    real*8 :: linear_momentum(1,3), linear_momentum_new(1,3)
    real*8 :: lpmag, lpmag_new
    ! angular momentum 
    real*8 :: j(1,3), j_new(1,3), r(ctrl%natom,3)
    real*8 :: com(3), jmag, jmag_new, summass
    ! decide if kick into cycle
    integer :: kick_into_cycle(ctrl%lpzpe_nah)

    if (printlevel>3) then
      write(u_log,*) '============================================================='
      write(u_log,*) '       Local Pair - Zero Point Energy Correction Scheme'
      write(u_log,*) '============================================================='
      write(u_log,*) 'S. Mukherjee, M. Barbatti. J. Chem. Theory Comput. 2022, 18, 4109-4116'
      write(u_log,*) 'Y. Shu, D. G. Truhlar, to be submitted'
    endif

    ! the big do loop over all AH bonds
    do ibond=1, ctrl%lpzpe_nah

      if (ctrl%lpzpe_base==0 .or. ctrl%lpzpe_base==1) then 
        kick_into_cycle(ibond)=1
      else if (ctrl%lpzpe_base==2) then 
        if (traj%lpzpe_iter_outcycle(ibond)>=2) then 
          kick_into_cycle(ibond)=1
        else
          kick_into_cycle(ibond)=0
        endif 
        traj%lpzpe_iter_outcycle(ibond)=traj%lpzpe_iter_outcycle(ibond)+1
      endif 

      if (kick_into_cycle(ibond)==1) then 

      if (ctrl%lpzpe_base==2 .and. traj%in_cycle(ibond)==0 .and. traj%lpzpe_iter_incycle(ibond)==1) then
        ! first iteration in a cummulative cycle for direct LP-ZPE, compute all properties
        ! compute the following properties on the fly 
        ! traj%lpzpe_ke_zpe_ah(ibond): from local mode ZPE
        ! traj%t_cycle(ibond): from local mode frequency
        atom1=ctrl%lpzpe_ah(ibond,1)
        atom2=ctrl%lpzpe_ah(ibond,2)
        call compute_dlp_zpe(ctrl%natom, traj%geom_ad, traj%geom_old_ad, traj%grad_ad, traj%grad_old_ad, &
          &traj%mass_a, atom1, atom2, traj%lpzpe_ke_zpe_ah(ibond), traj%t_cycle(ibond))
        traj%t_check(ibond)=traj%t_cycle(ibond) + ctrl%dtstep 
      endif

      if (ctrl%lpzpe_base==0 .or. ctrl%lpzpe_base==1) then
        tstart(ibond)=traj%lpzpe_endtime(ibond)+ctrl%dtstep
        tend(ibond)=tstart(ibond)+traj%t_cycle(ibond)
        !tstart(ibond)=traj%lpzpe_cycle(ibond)*traj%t_cycle(ibond)
        !tend(ibond)=traj%lpzpe_cycle(ibond)*traj%t_cycle(ibond) + traj%t_check(ibond)
      ! the time cycle varies. 
      else if (ctrl%lpzpe_base==2) then 
        tstart(ibond)=traj%lpzpe_endtime(ibond)+2*ctrl%dtstep
        tend(ibond)=tstart(ibond)+traj%t_cycle(ibond) 
      endif  

      ! within this cumulation period, computing average value of kinetic energy
      if (traj%microtime>=tstart(ibond) .and. traj%microtime<tend(ibond)) then
        if (traj%lpzpe_iter_incycle(ibond)==1) then
          !traj%lpzpe_starttime(ibond)=traj%microtime-ctrl%dtstep
          traj%lpzpe_starttime(ibond)=traj%microtime
          traj%in_cycle(ibond)=1
          if (printlevel>4) then
            write(u_log,*) '=====AH bond index:', ibond
            write(u_log,*) 'At the start of a cycle'
            write(u_log,'(a,F9.4,a)') 'local mode frequency:', 2*traj%lpzpe_ke_zpe_ah(ibond)*au2rcm, ' cm^-1'
            write(u_log,'(a,F9.4,a)') 'local mode zero-point energy:', traj%lpzpe_ke_zpe_ah(ibond)*au2eV, ' eV'
            write(u_log,'(a,F9.4,a)') 'local mode vibrational cycle:', traj%t_cycle(ibond)*au2fs, ' fs'
            write(u_log,'(a,F9.4,a)') 'Cycle start time:', tstart(ibond)*au2fs, ' fs'
            write(u_log,'(a,F9.4,a)') 'Cycle end time:', tend(ibond)*au2fs, ' fs'
            write(u_log,'(a,F9.4,a)') 'Cummulation of kinetic energy begins at:', traj%lpzpe_starttime(ibond)*au2fs, ' fs'
            write(u_log,*) 'Reset cumulative AH and BC bonds kinetic energy to zero'
          endif
          traj%lpzpe_ke_ah(ibond)=0.d0
        endif
        if (printlevel>4) then
          write(u_log,*) 'Within a LP-ZPE cycle, current cycle:', traj%lpzpe_cycle(ibond)
          write(u_log,*) 'current ieration in cycle:', traj%lpzpe_iter_incycle(ibond)
          write(u_log,*) 'cummulative AH bond kinetic energy of last step:', traj%lpzpe_ke_ah(ibond)
        endif
        ! cumulate ith AH bond kinetic energy
        atom1=ctrl%lpzpe_ah(ibond,1)
        atom2=ctrl%lpzpe_ah(ibond,2)
        call compute_parallel_velocity(ctrl%natom, traj%geom_ad, traj%veloc_ad, traj%mass_a, atom1, atom2, &
          &parallel_vec, parallel_veloc1, parallel_veloc2, ke_bond_ah)
        if (traj%lpzpe_iter_incycle(ibond)>1) then 
          traj%lpzpe_ke_ah(ibond)=traj%lpzpe_ke_ah(ibond)+ke_bond_ah*ctrl%dtstep
        endif
        if (printlevel>4) write(u_log,*) 'cummulative AH bond kinetic energy:', traj%lpzpe_ke_ah(ibond)
        traj%lpzpe_iter_incycle(ibond)=traj%lpzpe_iter_incycle(ibond)+1
      endif

      !cummulation is done, start computing averaged ZPE, and make ZPE correction
      if (traj%microtime>=tend(ibond) .and. traj%in_cycle(ibond)==1) then

        traj%lpzpe_endtime(ibond)=traj%microtime
        ! compute average kinetic energy
        traj%lpzpe_ke_ah(ibond)=traj%lpzpe_ke_ah(ibond)/(traj%lpzpe_endtime(ibond)-traj%lpzpe_starttime(ibond))
        ! reassign the base line as an averaged kinetic energy from first cycle
        ! Notice that this is original scheme used in LP-ZPE paper. 
        if (ctrl%lpzpe_base==0 .and. traj%lpzpe_cycle(ibond)==1) then
          if (printlevel>4) write(u_log,*) "Use averaged kinetic energy from first cycle as base line"
          traj%lpzpe_ke_zpe_ah(ibond)=traj%lpzpe_ke_ah(ibond)
        endif 
        ! print information. 
        if (printlevel>4) then
          write(u_log,*) 'At the end of a cycle'
          write(u_log,'(a,F9.4,a)') 'Cycle start time:', tstart(ibond)*au2fs, ' fs'
          write(u_log,'(a,F9.4,a)') 'Cycle end time:', tend(ibond)*au2fs, ' fs'
          write(u_log,'(a,F9.4,a)') 'Cummulation of kinetic energy  begins at:', traj%lpzpe_starttime(ibond)*au2fs, ' fs'
          write(u_log,'(a,F9.4,a)') 'Cummulation of kinetic energy  ends at:', traj%lpzpe_endtime(ibond)*au2fs, ' fs'
          write(u_log,'(a,F9.4)') 'Average AH bond kinetic energy', traj%lpzpe_ke_ah(ibond)
          call vec3write(ctrl%natom, traj%veloc_ad, u_log, ' velocity before LP-ZPE correction', 'F14.9')
          ! momentum
          do iatom=1,ctrl%natom
            momentum(iatom,:)=traj%mass_a(iatom)*traj%veloc_ad(iatom,:)
          enddo
          call vec3write(ctrl%natom, momentum, u_log, ' momentum before LP-ZPE correction', 'F14.9')
          ! center of mass momentum 
          linear_momentum(:,:)=0.d0
          do idir=1,3
            do iatom=1,ctrl%natom
              linear_momentum(1,idir)=linear_momentum(1,idir)+traj%mass_a(iatom)*traj%veloc_ad(iatom,idir)
            enddo
          enddo
          lpmag=0.d0
          do idir=1,3
            lpmag=lpmag+linear_momentum(1,idir)**2
          enddo
          lpmag=dsqrt(lpmag)
          call vec3write(1, linear_momentum, u_log, ' linear momentum before LP-ZPE correction', 'F14.9')
          write(u_log,'(A54,X,F14.9)') 'magnitude of linear momentum before LP-ZPE correction', lpmag
          ! angular momentum 
          j(:,:)=0.d0
          jmag=0.d0
          com(:)=0.d0
          summass=0.d0
          do iatom=1,ctrl%natom
            summass=summass+traj%mass_a(iatom)
            do idir=1,3
              com(idir)=com(idir)+traj%geom_ad(iatom,idir)*traj%mass_a(iatom)
            enddo
          enddo
          com=com/summass
          do iatom=1,ctrl%natom
            do idir=1,3
              r(iatom,idir)=traj%geom_ad(iatom,idir)-com(idir)
            enddo
          enddo
          do iatom=1,ctrl%natom
            j(1,1)=j(1,1)+(r(iatom,2)*momentum(iatom,3)-r(iatom,3)*momentum(iatom,2))
            j(1,2)=j(1,2)-(r(iatom,1)*momentum(iatom,3)-r(iatom,3)*momentum(iatom,1))
            j(1,3)=j(1,3)+(r(iatom,1)*momentum(iatom,2)-r(iatom,2)*momentum(iatom,1))
          enddo
          do idir=1,3
            jmag=jmag+j(1,idir)**2
          enddo
          jmag=dsqrt(jmag)
          call vec3write(1, j, u_log, ' angular momentum before LP-ZPE correction', 'F14.9')
          write(u_log,'(A54,X,F14.9)') 'magnitude of angular momentum before LP-ZPE correction', jmag
          ! kinetic energy
          Ekin=Calculate_ekin(ctrl%natom, traj%veloc_ad, traj%mass_a)
          write(u_log,'(A39,X,F14.9)') 'kinetic energy before LP-ZPE correction', Ekin
          ! bond energy 
          atom1=ctrl%lpzpe_ah(ibond,1)
          atom2=ctrl%lpzpe_ah(ibond,2)
          call compute_parallel_velocity(ctrl%natom, traj%geom_ad, traj%veloc_ad, traj%mass_a, atom1, atom2, &
            &parallel_vec, parallel_veloc1, parallel_veloc2, bond_Ekin_AH)
          write(u_log,'(A39,X,F14.9)') 'Kinetic energy of AH bond before LP-ZPE correction', bond_Ekin_AH
          write(u_log,*) '====== Start LP-ZPE correction ======'
        endif
        ! ===perform correction
        if (ctrl%lpzpe_scheme==0) then ! use the original scheme
          if (printlevel>4) write(u_log,*) 'Employ the original LP-ZPE scheme by Mukherjee and Barbatti'
          call LP_ZPE_MB(ctrl%natom, traj%geom_ad, traj%veloc_ad, traj%mass_a, &
            &traj%element_a, ctrl%lpzpe_nbc, ctrl%lpzpe_bc, ctrl%lpzpe_nah, ctrl%lpzpe_ah, &
            &ibond, traj%lpzpe_ke_zpe_ah(ibond), traj%lpzpe_ke_ah(ibond), ctrl%ke_threshold)
        else if (ctrl%lpzpe_scheme==1) then 
          if (printlevel>4) write(u_log,*) 'Employ the improved LP-ZPE scheme by Shu and Truhlar'
          call LP_ZPE_ST(ctrl%natom, traj%geom_ad, traj%veloc_ad, traj%mass_a, &
            &traj%element_a, ctrl%lpzpe_nbc, ctrl%lpzpe_bc, ctrl%lpzpe_nah, ctrl%lpzpe_ah, &
            &ibond, traj%lpzpe_ke_zpe_ah(ibond), traj%lpzpe_ke_ah(ibond), ctrl%ke_threshold, &
            &ctrl%ilpzpe_f1, ctrl%ilpzpe_f2, ctrl%ilpzpe_f3, ctrl%ilpzpe_f4)
        endif 
        ! print information
        if (printlevel>4) then 
          write(u_log,*) '====== LP-ZPE correction finished ======'
          call vec3write(ctrl%natom, traj%veloc_ad, u_log, ' velocity after LP-ZPE correction', 'F14.9')
          ! momentum 
          do iatom=1,ctrl%natom
            momentum_new(iatom,:)=traj%mass_a(iatom)*traj%veloc_ad(iatom,:)
          enddo
          call vec3write(ctrl%natom, momentum_new, u_log, ' momentum after LP-ZPE correction', 'F14.9')
          ! center of mass momentum 
          linear_momentum_new(:,:)=0.d0
          do idir=1,3
            do iatom=1,ctrl%natom
              linear_momentum_new(1,idir)=linear_momentum_new(1,idir)+traj%mass_a(iatom)*traj%veloc_ad(iatom,idir)
            enddo
          enddo
          lpmag_new=0.d0
          do idir=1,3
            lpmag_new=lpmag_new+linear_momentum(1,idir)**2
          enddo
          lpmag_new=dsqrt(lpmag_new)
          call vec3write(1, linear_momentum_new, u_log, ' linear momentum after LP-ZPE correction', 'F14.9')
          write(u_log,'(A54,X,F14.9)') 'magnitude of linear momentum after LP-ZPE correction', lpmag_new
          ! angular momentum 
          j_new(:,:)=0.d0
          jmag_new=0.d0
          do iatom=1,ctrl%natom
            j_new(1,1)=j_new(1,1)+(r(iatom,2)*momentum_new(iatom,3)-r(iatom,3)*momentum_new(iatom,2))
            j_new(1,2)=j_new(1,2)-(r(iatom,1)*momentum_new(iatom,3)-r(iatom,3)*momentum_new(iatom,1))
            j_new(1,3)=j_new(1,3)+(r(iatom,1)*momentum_new(iatom,2)-r(iatom,2)*momentum_new(iatom,1))
          enddo
          do idir=1,3
            jmag_new=jmag_new+j_new(1,idir)**2
          enddo
          jmag_new=dsqrt(jmag_new)
          call vec3write(1, j_new, u_log, ' angular momentum after LP-ZPE correction', 'F14.9')
          write(u_log,'(A54,X,F14.9)') 'magnitude of angular momentum after LP-ZPE correction', jmag_new
          ! kinetic energy
          Ekin_new=Calculate_ekin(ctrl%natom, traj%veloc_ad, traj%mass_a)
          write(u_log,'(A39,X,F14.9)') 'kinetic energy after LP-ZPE correction', Ekin_new
          ! bond energy 
          atom1=ctrl%lpzpe_ah(ibond,1)
          atom2=ctrl%lpzpe_ah(ibond,2)
          call compute_parallel_velocity(ctrl%natom, traj%geom_ad, traj%veloc_ad, traj%mass_a, atom1, atom2, &
            &parallel_vec, parallel_veloc1, parallel_veloc2, bond_Ekin_AH_new)
          write(u_log,'(A39,X,F14.9)') 'Kinetic energy of AH bond after LP-ZPE correction', bond_Ekin_AH_new
        endif 
        ! kick trajectory out of a cycle
        traj%in_cycle(ibond)=0
        traj%lpzpe_ke_ah(ibond)=0.d0
        traj%lpzpe_iter_incycle(ibond)=1
        traj%lpzpe_iter_outcycle(ibond)=1
        traj%lpzpe_cycle(ibond)=traj%lpzpe_cycle(ibond)+1

      endif !if (traj%microtime>=tend .and. traj%in_cycle==1) then

      endif ! if (kick_into_cycle(ibond)==1) then

    enddo ! the big do loop over AH bonds.

  endsubroutine 

! ==========================================================
!> subroutine to compute parallel velocity
  subroutine compute_parallel_velocity(n, g, v, m, atom1, atom2, parallel_vec, &
    &parallel_veloc1, parallel_veloc2, ke)
    implicit none
    integer, intent(in) :: n
    real*8, intent(in) :: g(n,3), v(n,3), m(n)
    integer, intent(in) :: atom1, atom2
    real*8, intent(inout) :: parallel_vec(3)
    real*8, intent(inout) :: parallel_veloc1, parallel_veloc2
    real*8, intent(inout) :: ke
    real*8 :: r_m
    real*8 :: norm_parallel
    real*8 :: pair_veloc
    integer :: idir

    parallel_veloc1=0.d0
    parallel_veloc2=0.d0

    r_m=m(atom1)*m(atom2)/(m(atom1)+m(atom2))

    parallel_vec(:)=g(atom1,:)-g(atom2,:)

    norm_parallel=0.d0
    do idir=1,3
      norm_parallel=norm_parallel+parallel_vec(idir)*parallel_vec(idir)
    enddo
    norm_parallel=sqrt(norm_parallel)
    parallel_vec=parallel_vec/norm_parallel

    do idir=1,3
      parallel_veloc1=parallel_veloc1+v(atom1,idir)*parallel_vec(idir)
      parallel_veloc2=parallel_veloc2+v(atom2,idir)*parallel_vec(idir)
    enddo

    pair_veloc=parallel_veloc1-parallel_veloc2

    ke=0.5*r_m*pair_veloc*pair_veloc

  endsubroutine

! ==========================================================
!> original LP-ZPE scheme by Mukherjee and Barbatti
!> LP_ZPE_MB(ctrl%natom, traj%geom_ad, traj%veloc_ad, traj%mass_a, &
!>   &traj%element_a, ctrl%lpzpe_nbc, ctrl%lpzpe_bc, ctrl%lpzpe_nah, ctrl%lpzpe_ah, &
!>   &ibond, traj%lpzpe_ke_zpe_ah(ibond), traj%lpzpe_ke_ah(ibond), ctrl%ke_threshold)
  subroutine LP_ZPE_MB(n, g, v, m, element_a, nbc, list_bc, nah, list_ah, &
    &ibond, ah_zpe, ke_ah, ke_threshold)
    use definitions, only: u_log, printlevel
    use matrix

    implicit none
    integer, intent(in) :: n
    real*8, intent(in) :: g(n,3), m(n)
    character*2, intent(in) :: element_a(n)
    integer, intent(in) :: nbc, list_bc(nbc,2)
    integer, intent(in) :: nah, list_ah(nah,2)
    integer, intent(in) :: ibond
    real*8, intent(in) :: ah_zpe
    real*8, intent(in) :: ke_ah, ke_threshold
    real*8, intent(inout) :: v(n,3)

    integer :: do_correction ! 0=not do, 1=do
    real*8 :: e_leak
    integer :: scheme
    real*8 :: r_m_ah
    real*8 :: parallel_vec(nbc,3)
    real*8 :: vB(nbc), vC(nbc)
    integer :: atom1, atom2
    ! for computing BC corrections
    integer :: atom3, atom4
    real*8 :: ke(nbc), sum_ke, sum_ke_new
    integer :: bc_valid(nbc)
    integer :: rescale
    real*8 :: factor(nbc), deltaE(nbc)
    integer :: iter, maxiter, remaining_nbc

    integer :: jbond, iatom, jatom, katom, idir, ipair

    atom1=list_ah(ibond,1)
    atom2=list_ah(ibond,2)
    e_leak=ah_zpe-ke_ah
    r_m_ah=m(atom1)*m(atom2)/(m(atom1)+m(atom2))

    if (e_leak>ke_threshold) then
      if (printlevel>4) then
        write(u_log,'(A25,I3,X,A2,X,A2,X,A8)') ' ZPE leaking for AH bond ', &
          &ibond, element_a(atom1), element_a(atom2), 'detected'
        write(u_log,'(A35,X,F14.9,X,A10,F14.9)') '    current AH bond kinetic energy=', &
          &ke_ah, 'should be=', ah_zpe
        write(u_log,'(A19,F14.9,X,A10,F14.9)') '    leaking energy=', e_leak, &
          &'threshold=', ke_threshold
      endif
      do_correction=1
      ! correction on BC bonds and check if we have enough kinetic energy
      if (printlevel>4) write(u_log,*) ' --Start correcting BC bonds--'
      ! ==start computing corrections on each BC bond
      maxiter=nbc
      iter=1
      remaining_nbc=nbc
      ! initialize all BC bonds to be valid in computing factor
      bc_valid(:)=1
      do jbond=1,nbc
        atom3=list_bc(jbond,1)
        atom4=list_bc(jbond,2)
        call compute_parallel_velocity(n, g, v, m, atom3, atom4, parallel_vec(jbond,:), &
          &vB(jbond), vC(jbond), ke(jbond))
      enddo
      rescale=1 ! so that we start with the do while loop
      do while (iter.le.maxiter .and. rescale.eq.1)
        rescale=0
        sum_ke=0.d0
        factor(:)=0.d0
        do jbond=1,nbc
          if (bc_valid(jbond)==1) then
            sum_ke=sum_ke+ke(jbond)
          endif
        enddo
        do jbond=1,nbc
          if (bc_valid(jbond)==1) then
            factor(jbond)=ke(jbond)/sum_ke
            deltaE(jbond)=(vB(jbond)-vC(jbond))**2-(2*e_leak*factor(jbond)/r_m_ah)
            if (deltaE(jbond)<0) then
              bc_valid(jbond)=0
              rescale=1
              remaining_nbc=remaining_nbc-1
            endif
          endif
        enddo
        iter=iter+1
      enddo
      if (remaining_nbc==0) then
        do_correction=0
      else if (remaining_nbc>0) then
        ! correct BC velocity
        scheme=2
        do jbond=1,nbc
          atom3=list_bc(jbond,1)
          atom4=list_bc(jbond,2)
          if (bc_valid(jbond)==1) then
            call correct_bond_velocity(n, g, v, m, element_a, atom3, atom4, &
              &e_leak, factor(jbond), r_m_ah, scheme)
          endif
        enddo
      endif !if (remaining_nbc==0) then
      ! ==end computing corrections on each BC bond
      if (do_correction==1) then ! we have enough kinetic energy form BC bonds, use that to correct AH velocity
        ! correct AH velocity
        scheme=1
        if (printlevel>4) write(u_log,*) ' --Start correcting AH bonds--'
        call correct_bond_velocity(n, g, v, m, element_a, atom1, atom2, &
          &e_leak, 1.d0, r_m_ah, scheme)
      else if (do_correction==0) then
        if (printlevel>4) write(u_log,*) ' Skip correction, becasue BC bonds do not have enough kinetic energy'
      endif
    endif ! if (e_leak>ke_threshold) then

  endsubroutine

! ==========================================================
!> subroutine to correct bond velocity
!> correct_bond_velocity(ctrl%natom, traj%geom_ad, traj%veloc_ad, traj%mass_a, &
!>   &atom1, atom2, eleak, factor, scheme)
  subroutine correct_bond_velocity(n, g, v, m, element_a, atom1, atom2, e, factor, r_m_ah, scheme)
    implicit none
    integer, intent(in) :: n
    real*8, intent(in) :: g(n,3), m(n)
    character*2, intent(in) :: element_a(n)
    integer, intent(in) :: atom1, atom2
    real*8, intent(in) :: e, factor, r_m_ah
    integer, intent(in) :: scheme
    real*8, intent(inout) :: v(n,3)

    real*8 :: parallel_vec(3)
    real*8 :: parallel_veloc1, parallel_veloc2
    real*8 :: pair_veloc, r_m, ke_old, ke_new
    real*8 :: deltaV
    integer :: idir

    call compute_parallel_velocity(n, g, v, m, atom1, atom2, parallel_vec, &
      &parallel_veloc1, parallel_veloc2, ke_old)

    r_m=m(atom1)*m(atom2)/(m(atom1)+m(atom2))
    pair_veloc=parallel_veloc1-parallel_veloc2

    if (scheme==1) then
      deltaV=sqrt(pair_veloc**2+(2*e*factor/r_m_ah))-pair_veloc
    else if (scheme==2) then
      deltaV=sqrt(pair_veloc**2-(2*e*factor/r_m))-pair_veloc
    endif

    do idir=1,3
      v(atom1,idir)=v(atom1,idir) + r_m/m(atom1)*deltaV*(parallel_vec(idir))
      v(atom2,idir)=v(atom2,idir) - r_m/m(atom2)*deltaV*(parallel_vec(idir))
    enddo

  endsubroutine


! ==========================================================
!> improved LP-ZPE scheme by Shu and Truhlar
!> LP_ZPE_ST(ctrl%natom, traj%geom_ad, traj%veloc_ad, traj%mass_a, &
!>  &traj%element_a, ctrl%lpzpe_nbc, ctrl%lpzpe_bc, ctrl%lpzpe_nah, ctrl%lpzpe_ah, &
!>  &traj%lpzpe_ke_zpe_ah(ibond), traj%lpzpe_ke_ah(ibond), ctrl%ke_threshold, &
!>  &ctrl%ilpzpe_f1, ctrl%ilpzpe_f2, ctrl%ilpzpe_f3, ctrl%ilpzpe_f4
  subroutine LP_ZPE_ST(n, g, v, m, element_a, nbc, list_bc, &
    &nah, list_ah, ibond, ah_zpe, ke_ah, ke_threshold, f1, f2, f3, f4)
    use definitions, only: u_log, printlevel
    use matrix

    implicit none
    integer, intent(in) :: n
    real*8, intent(in) :: g(n,3), m(n)
    character*2, intent(in) :: element_a(n)
    integer, intent(in) :: nbc, list_bc(nbc,2)
    integer, intent(in) :: nah, list_ah(nah,2)
    integer, intent(in) :: ibond
    real*8, intent(in) :: ah_zpe, ke_ah, ke_threshold
    real*8, intent(in) :: f1, f2, f3, f4
    real*8, intent(inout) :: v(n,3)

    real*8 :: e_leak
    real*8 :: parallel_vec(3)
    real*8 :: parallel_veloc1
    real*8 :: parallel_veloc2
    ! used for saving all energy correction 
    real*8 :: bond_Ekin_AH
    ! angular momentum 
    real*8 :: r(n,3), com(3), summass
    ! for computing bvec (a single number for single AH bond)
    real*8 :: dveloc, dEvec, bvec, rm
    ! for optimization 
    real*8 :: opt_value, grad_opt_value(n,3)
    ! for check 
    real*8 :: bvec_computed, Ecorrection_computed
    ! final correction
    real*8 :: veloc_correction(n,3)

    integer :: atom1, atom2
    integer :: jbond, iatom, jatom, katom, idir
 
    e_leak=ah_zpe-ke_ah
   
    if (e_leak>ke_threshold) then ! Start performing correction. 

      atom1=list_ah(ibond,1)
      atom2=list_ah(ibond,2)
      call compute_parallel_velocity(n, g, v, m, atom1, atom2, &
        &parallel_vec, parallel_veloc1, parallel_veloc2, bond_Ekin_AH)
      rm=m(atom1)*m(atom2)/(m(atom1)+m(atom2))
      dveloc=parallel_veloc1-parallel_veloc2
      dEvec=dveloc**2+2*e_leak/rm
      bvec=(-dveloc+sqrt(dEvec))
      ! conservation of angular momentum 
      com(:)=0.d0
      summass=0.d0
      do iatom=1,n
        summass=summass+m(iatom)
        do idir=1,3
          com(idir)=com(idir)+g(iatom,idir)*m(iatom)
        enddo
      enddo
      com=com/summass
      do iatom=1,n
        do idir=1,3
          r(iatom,idir)=g(iatom,idir)-com(idir)
        enddo
      enddo
      ! initialize correction on velocity 
      veloc_correction(:,:)=0.d0
      call random_number(veloc_correction)
      veloc_correction=1000*veloc_correction
      ! optimization of velocity correction term
      call optimize_veloc_correction(n, r, v, m, veloc_correction, &
        &nbc, list_bc, nah, list_ah, ibond, parallel_vec, bvec, f1, f2, f3, f4)
      ! now veloc_correction is optimized, check if it is positive correction
      bvec_computed=0.d0
      Ecorrection_computed=0.d0
      bvec_computed=0.d0
      do idir=1,3
        bvec_computed=bvec_computed+&
          &(veloc_correction(atom1,idir)-veloc_correction(atom2,idir))*parallel_vec(idir)
      enddo 
      Ecorrection_computed=rm*(parallel_veloc1-parallel_veloc2)*bvec_computed+0.5*rm*bvec_computed**2
      if (Ecorrection_computed<0.d0) then 
        if (printlevel>4) write(u_log,*) "total correction to AH bonds are negative, set correction to zero"
        veloc_correction(:,:)=0.d0
      endif  
      if (printlevel>4) then
        write(u_log,'(A7,X,I3,X,A18,F14.9,A18,F14.9)') 'AH bond', ibond, &
          &'corrected amount: ', Ecorrection_computed, 'required amount: ', e_leak
      endif
      ! addition of correction.
      do iatom=1,n
        do idir=1,3
          v(iatom,idir)=v(iatom,idir)+veloc_correction(iatom,idir)
        enddo
      enddo

    else ! no leak, exit if, and no correction is performed. 
      if (printlevel>4) then
        write(u_log,*) 'No leakage found, skip correction'
      endif
    endif ! if (e_leak>ke_threshold) then

  endsubroutine

! ==========================================================
!> subroutine to optimize veloc_correction using L-BFGS
!> subroutine optimize_veloc_correction(n, r, v, m, veloc_correction, &
!>   &nbc, list_bc, nah, list_ah, ibond, parallel_vec, bvec, f1, f2, f3, f4)
  subroutine optimize_veloc_correction(n, r, v, m, veloc_correction, &
    &nbc, list_bc, nah, list_ah, ibond, parallel_vec, bvec, f1, f2, f3, f4)
    use definitions, only: u_log, printlevel
    use matrix
    use optimizer

    implicit none
    integer, intent(in) :: n
    real*8, intent(in) :: r(n,3) !relative geometri w.r.t. center of mass
    real*8, intent(in) :: v(n,3)
    real*8, intent(in) :: m(n)
    real*8, intent(inout) :: veloc_correction(n,3)
    integer, intent(in) :: nbc, list_bc(nbc,2)
    integer, intent(in) :: nah, list_ah(nah,2)
    integer, intent(in) :: ibond
    real*8, intent(in) :: parallel_vec(3)
    real*8, intent(in) :: bvec
    real*8, intent(in) :: f1, f2, f3, f4
    real*8 :: opt_value, grad_opt_value(n,3)

    integer :: opt_n, opt_m, iprint
    real*8, parameter :: factr  = 1.0e+1_wp, pgtol  = 1.0e-5_wp
    character(len=60) :: task, csave
    logical :: lsave(4)
    integer :: isave(44)
    real*8 :: f, dsave(29)
    integer, allocatable :: nbd(:), iwa(:)
    real*8, allocatable :: x(:), l(:), u(:), g(:), wa(:)
    real*8 :: veloc_correction_tmp(n,3)

    integer :: i, j, k

    opt_n=3*n
    opt_m=5
    if (printlevel>3) then
      iprint=1
    else
      iprint=0
    endif
    allocate (nbd(opt_n), x(opt_n), l(opt_n), u(opt_n), g(opt_n))
    allocate (iwa(3*opt_n))
    allocate (wa(2*opt_m*opt_n + 5*opt_n + 11*opt_m*opt_m + 8*opt_m) )

    do i=1,opt_n
      nbd(i)=0
    enddo

    ! convert veloc_correction to x 
    do i=1,n
      do j=1,3
        k=(i-1)*3+j
        x(k)=veloc_correction(i,j)
      enddo
    enddo

    veloc_correction_tmp(:,:)=veloc_correction(:,:)

    task='START'
    do while(task(1:2)=='FG'.or.task=='NEW_X'.or.task=='START')

      call setulb(opt_n, opt_m, x, l, u, nbd, f, g, factr, pgtol, &
        &wa, iwa, task, iprint, csave, lsave, isave, dsave)
      if (task(1:2) == 'FG') then
        do i=1,n
          do j=1,3
            k=(i-1)*3+j
            veloc_correction_tmp(i,j)=x(k)
          enddo
        enddo
        call objective_function(n, r, v, m, veloc_correction_tmp, &
          &nbc, list_bc, nah, list_ah, ibond, parallel_vec, bvec, f1, f2, f3, f4, &
          &opt_value, grad_opt_value)
        f=opt_value
        do i=1,n
          do j=1,3
            k=(i-1)*3+j
            g(k)=grad_opt_value(i,j)
          enddo
        enddo
      endif !if (task(1:2) == 'FG') then
    enddo

    veloc_correction(:,:)=veloc_correction_tmp(:,:)

  endsubroutine

! ==========================================================
! objective function used in optimization of veloc_correction
! call objective_function(n, v, m, veloc_correction, &
!   &nbc, list_bc, nah, list_ah, ibond parallel_vec, bvec, f1, f2, f3, f4, &
!   &opt_value, grad_opt_value)
  subroutine objective_function(n, r, v, m, veloc_correction, nbc, &
    &list_bc, nah, list_ah, ibond, parallel_vec, bvec, f1, f2, f3, f4, &
    &opt_value, grad_opt_value)
    use definitions, only: u_log, printlevel
    use matrix
    implicit none
    integer, intent(in) :: n
    real*8, intent(in) :: r(n,3) !relative geometri w.r.t. center of mass
    real*8, intent(in) :: v(n,3)
    real*8, intent(in) :: m(n)
    real*8, intent(in) :: veloc_correction(n,3)
    integer, intent(in) :: nbc, list_bc(nbc,2)
    integer, intent(in) :: nah, list_ah(nah,2)
    integer, intent(in) :: ibond
    real*8, intent(in) :: parallel_vec(3)
    real*8, intent(in) :: bvec
    real*8, intent(in) :: f1, f2, f3, f4
    real*8, intent(inout) :: opt_value, grad_opt_value(n,3)

    integer :: atom1, atom2
    real*8 :: sum_ke, cdotu
    real*8 :: deltap(3)
    integer :: iatom, idir
    ! angular momentum 
    real*8 :: deltaj(3), jmag

    opt_value=0.d0
    grad_opt_value(:,:)=0.d0

    atom1=list_ah(ibond,1)
    atom2=list_ah(ibond,2)
    cdotu=0.d0
    do idir=1,3
      cdotu=cdotu+&
        &(veloc_correction(atom1,idir)-veloc_correction(atom2,idir))*parallel_vec(idir)
    enddo
    opt_value=opt_value+f1*(cdotu-bvec)**2
    do idir=1,3
      grad_opt_value(atom1,idir)=grad_opt_value(atom1,idir)+&
        &2*f1*(cdotu-bvec)*parallel_vec(idir)
      grad_opt_value(atom2,idir)=grad_opt_value(atom2,idir)-&
        &2*f1*(cdotu-bvec)*parallel_vec(idir)
    enddo
    ! conservation of total kinetic energy
    sum_ke=0.d0
    do iatom=1,n
      do idir=1,3
        sum_ke=sum_ke+m(iatom)*v(iatom,idir)*veloc_correction(iatom,idir)+&
          &0.5*m(iatom)*veloc_correction(iatom,idir)*veloc_correction(iatom,idir)
      enddo
    enddo
    opt_value=opt_value+f2*sum_ke**2
    do iatom=1,n
      do idir=1,3
        grad_opt_value(iatom,idir)=grad_opt_value(iatom,idir)+&
          &2*f2*sum_ke*(m(iatom)*v(iatom,idir)+m(iatom)*veloc_correction(iatom,idir))
      enddo
    enddo
    ! conservation of linear momentum 
    deltap(:)=0.d0
    do idir=1,3
      do iatom=1,n
        deltap(idir)=deltap(idir)+m(iatom)*veloc_correction(iatom,idir)
      enddo
    enddo
    do idir=1,3
      opt_value=opt_value+f3*deltap(idir)**2
    enddo
    do idir=1,3
      do iatom=1,n
        grad_opt_value(iatom,idir)=grad_opt_value(iatom,idir)+&
          &2*f3*deltap(idir)*m(iatom)
      enddo
    enddo
    ! conservation of angular momentum 
    deltaj(:)=0.d0
    do iatom=1,n
      deltaj(1)=deltaj(1)+r(iatom,2)*m(iatom)*veloc_correction(iatom,3)-&
        &r(iatom,3)*m(iatom)*veloc_correction(iatom,2)
      deltaj(2)=deltaj(2)+r(iatom,1)*m(iatom)*veloc_correction(iatom,3)-&
        &r(iatom,3)*m(iatom)*veloc_correction(iatom,1)
      deltaj(3)=deltaj(3)+r(iatom,1)*m(iatom)*veloc_correction(iatom,2)-&
        &r(iatom,2)*m(iatom)*veloc_correction(iatom,1)
    enddo
    do idir=1,3
      opt_value=opt_value+f4*deltaj(idir)**2
    enddo
    do iatom=1,n
      grad_opt_value(iatom,3)=grad_opt_value(iatom,3)+&
        &2*f4*deltaj(1)*r(iatom,2)*m(iatom)
      grad_opt_value(iatom,2)=grad_opt_value(iatom,2)-&
        &2*f4*deltaj(1)*r(iatom,3)*m(iatom)
      grad_opt_value(iatom,3)=grad_opt_value(iatom,3)+&
        &2*f4*deltaj(2)*r(iatom,1)*m(iatom)
      grad_opt_value(iatom,1)=grad_opt_value(iatom,1)-&
        &2*f4*deltaj(2)*r(iatom,3)*m(iatom)
      grad_opt_value(iatom,2)=grad_opt_value(iatom,2)+&
        &2*f4*deltaj(3)*r(iatom,1)*m(iatom)
      grad_opt_value(iatom,1)=grad_opt_value(iatom,1)-&
        &2*f4*deltaj(3)*r(iatom,2)*m(iatom)
    enddo

  endsubroutine

! ==========================================================
!> compute ZPE of AH bond on the fly
  subroutine compute_dlp_zpe(n, g, g_old, grad, grad_old, m, atom1, atom2, zpe, t_cycle)
    use definitions, only: au2rcm, au2fs

    implicit none
    integer, intent(in) :: n
    real*8, intent(in) :: g(n,3), g_old(n,3), grad(n,3), grad_old(n,3), m(n)
    integer, intent(in) :: atom1, atom2
    real*8, intent(inout) :: zpe, t_cycle

    real*8 :: parallel_vec(3), norm_parallel
    real*8 :: parallel_vec_old(3), norm_parallel_old
    real*8 :: grad_parallel, grad_parallel_old
    real*8 :: grad_parallel1, grad_parallel2
    real*8 :: grad_parallel_old1, grad_parallel_old2
    real*8 :: rm, dd2, freq

    integer :: idir

    ! parallel vector along a AH bond
    ! norm_parallel = d_AH, norm_parallel_old = d_AH_old
    parallel_vec(:)=g(atom1,:)-g(atom2,:)
    parallel_vec_old(:)=g_old(atom1,:)-g_old(atom2,:)
    norm_parallel=0.d0
    norm_parallel_old=0.d0
    do idir=1,3
      norm_parallel=norm_parallel+parallel_vec(idir)*parallel_vec(idir)
      norm_parallel_old=norm_parallel_old+parallel_vec_old(idir)*parallel_vec_old(idir)
    enddo
    norm_parallel=sqrt(norm_parallel)
    norm_parallel_old=sqrt(norm_parallel_old)
    parallel_vec=parallel_vec/norm_parallel
    parallel_vec_old=parallel_vec_old/norm_parallel_old

    ! grad_paralel = dE/d(d_AH); grad_paralel_old = dE_old/d(d_AH_old)
    grad_parallel=0.d0
    grad_parallel_old=0.d0
    do idir=1,3
      grad_parallel=grad_parallel+(grad(atom1,idir)-grad(atom2,idir))*parallel_vec(idir)
      grad_parallel_old=grad_parallel_old+(grad_old(atom1,idir)-grad_old(atom2,idir))*parallel_vec_old(idir)
    enddo

    rm=m(atom1)*m(atom2)/(m(atom1)+m(atom2))
    dd2=(grad_parallel-grad_parallel_old)/(norm_parallel-norm_parallel_old)
    freq=sqrt(dd2/rm)
    ! notice that this is 0.25 instead of 0.5 becasue this is not real ZPE
    ! inestad, we want half of ZPE. 
    zpe=0.25d0*freq
    t_cycle=1/(freq*au2rcm*2.997925)*10**5/au2fs

  endsubroutine

! ===========================================================
! End LP-ZPE algorithm
! ===========================================================

endmodule
