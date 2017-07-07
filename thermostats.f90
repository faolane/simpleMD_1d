! *************************************************
! simpleMD 1D -  Module file
! *************************************************
! Contains subroutines and functions for the
! different thermostats
! *************************************************
! Copyright (C) 2017 Fabien Brieuc under the terms
! of the GNU General Public License.
! *************************************************

module thermostats
   use param
   implicit none
   private
   public thermostating, langevin, qtb

contains

   subroutine thermostating(mode)
      ! initialize / terminate thermostat (e.g. generate random force for qtb,
      ! open, close and suppress the random force file)
      use param
      implicit none

      integer :: j,k

      character(len=*), intent(in) :: mode

      if (mode == 'initialize') then
         select case(therm)
            case('langevin')
               allocate(R(nrep))
               std = dsqrt(2.D0 * gam * M * KB * T / dt)
            case('qtb')
               allocate(R(nrep)) ! store random forces
               allocate(Rk(nrep)) ! store random forces on NM
               allocate(C(nrep,0:nrep-1)) ! Transformation matrix NM -> direct
               ! Init transformation matrix
               do k = 0, nrep - 1
                  do j = 1, nrep
                     if (k == 0) C(j,k) = dsqrt(ONREP)
                     if (2*k > 0 .and. 2*k < nrep) then
                        C(j,k) = dsqrt(2.d0*ONREP)*dcos(TWOPI*dfloat(k*j)*ONREP)
                     endif
                     if (2*k == nrep) C(j,k) = dsqrt(ONREP) * (-1.d0)**j
                     if (2*k > nrep) then
                        C(j,k) = dsqrt(2.d0*ONREP)*dsin(TWOPI*dfloat(k*j)*ONREP)
                     endif
                  enddo
               enddo
               ! Random force generation : Standard Dammak et al
               if (noiseGenMode == 1) then
                  ! new random force at every steps
                  mQTB = 1
                  std = dsqrt(2.D0 * gam * M / (mQTB * dt))
                  call generQTBRandF1('qtbRandF.tmp')
                  open(randfFileUnit, file ='qtbRandF.tmp', form='unformatted')
               ! Random force generation : Modified Dammak et al
               else if (noiseGenMode == 2) then
                  ! new random force at every mQTB = 1/(2*fcut*dt) steps
                  mQTB = nint(1.d0 / (2.d0 * fcut * dt))
                  if (mQTB == 0) mQTB = 1
                  print*, 'new random forces every', mQTB, 'steps'
                  std = dsqrt(2.D0 * gam * M / (mQTB * dt))
                  call generQTBRandF2('qtbRandF.tmp')
                  open(randfFileUnit, file ='qtbRandF.tmp', form='unformatted')
               ! Random force generation : Barrat and Rodney
               else if (noiseGenMode == 3) then
                  ! new random force at every mQTB = 1/(2*fcut*dt) steps
                  mQTB = nint(1.d0 / (2.d0 * fcut * dt))
                  if (mQTB == 0) mQTB = 1
                  print*, 'new random forces every', mQTB, 'steps'
                  std = dsqrt(2.D0 * gam * M / (mQTB * dt))
                  ! ensure that nQTB is even
                  if (mod(nQTB,2) /= 0) nQTB = nQTB + 1
                  allocate(HQTB(0:nQTB-1)) !"filter" in the time domain
                  allocate(rQTB(nrep,0:nQTB-1)) !random numbers
                  call generQTBRandF3()

               endif
         end select
      else if (mode == 'terminate') then
         select case(therm)
            case('langevin')
               deallocate(R)
            case('qtb')
               deallocate(C)
               deallocate(R)
               deallocate(Rk)
               if (noiseGenMode == 1 .or. noiseGenMode == 2) then
                  close(randfFileUnit)
                  call system('rm -f qtbRandF.tmp')
               else if (noiseGenMode == 3) then
                  deallocate(HQTB)
                  deallocate(rQTB)
               endif
         end select
      endif

   end subroutine thermostating

   real(dp) function langevin(v,irep)
      ! compute the forces associated with the standard (ie. classical)
      ! Langevin thermostat
      ! ******************************************
      ! The random fore is a gaussian white noise:
      ! <R(t)> = 0
      ! <R(t)R(t+tau)> = 2 * m * gam * kT * delta(tau)
      use param
      implicit none

      real(dp), intent(in) :: v
      integer, intent(in)  :: irep ! replica nb

      R(irep) = gaussianRand(std, 0.d0)
      langevin = R(irep) - m * gam * v

   end function langevin

   real(dp) function qtb(v, irep, n)
      ! compute forces associated with the QTB thermostat (ie. quantum
      ! Langevin thermostat)
      use param
      implicit none

      real(dp), intent(in) :: v ! velocity
      integer, intent(in)  :: n ! step number
      integer, intent(in)  :: irep ! replica nb
      integer :: i, j, k
      real :: R_sngl

      !new QTB random force every mQTB step
      if (mod(n, mQTB) == 0 .and. irep == 1) then
         ! read / generates forces for all normal modes
         if (noiseGenMode == 1 .or. noiseGenMode == 2) then
            do j = 1, nrep
               read(randfFileUnit) R_sngl
               Rk(j) = R_sngl
            enddo
         else if (noiseGenMode == 3) then
            do j = 1, nrep
               ! compute random force through the convolution product
               Rk(j) = 0.d0
               do i = 0, nQTB-1
                  Rk(j) = Rk(j) + HQTB(i) * rQTB(j, nQTB - i - 1)
               enddo
               ! update random numbers
               do i = 0, nQTB-2
                  rQTB(j, i) = rQTB(j, i+1)
               enddo
               rQTB(j, nQTB-1) = gaussianRand(1.d0,0.d0)
            enddo
            if (piqtbMode == 1 .and. nrep > 1) then
               ! centroid is classical
               Rk(1) = gaussianRand(dsqrt(KB*T),0.d0)
            endif
         endif
         ! Transform forces from NM to replicas
         do j = 1, nrep
            R(j) = 0.d0
            do k = 1, nrep
               R(j) = R(j) + C(j,k-1) * Rk(k)
            enddo
         enddo
      endif

      qtb = std * R(irep) - m * gam * v

   end function qtb

   subroutine generQTBRandF1(filename)
      ! generates forces associated with the QTB thermostat (ie. quantum
      ! Langevin thermostat)
      ! *********************************************
      ! if nrep = 1 : standard QTB thermostat
      ! *********************************************
      ! The random force is a stationnary gaussian colored noise with
      ! a Power Spectral Density (PSD) given by:
      ! 2 * m * gamma * theta(w, T) and
      ! theta(w, T) = hbar * w * ( 0.5 + 1 / (exp( hbar * w / kT ) - 1 ) )
      ! so that:
      ! <R(t)> = 0
      ! <R(t)R(t+tau) = FT(2 * m * gamma * theta(w,T)) (Wiener-Khinchin)
      ! see Dammak et al, PRL, 103, 190601 (2009) for more details
      ! **********************************************
      ! if nrep > 1 : QTB-PIMD thermostat
      ! **********************************************
      ! The random force is also a stationnary gaussian colored noise but with
      ! a modified Power Spectral Density (PSD) given by:
      ! 2 * m * gamma * kappa(w, T) and
      ! kappa(w, T) = P * kT * f_P(hw/2kT) and f_P solution of
      ! eq. (13) for f_P^(0) and solution of eq. (16) for f_P^(1)
      ! so that:
      ! <R(t)> = 0
      ! <R(t)R(t+tau) = FT(2 * m * gamma * kappa(w,T)) (Wiener-Khinchin)
      ! see Brieuc, Dammak and Hayoun, JCTC, 12, 1351-1359 (2016)
      ! **********************************************
      ! the random force is generated here using the method proposed by Dammak
      ! and coworkers and described in details in the appendix of
      ! Brieuc et al., JCTC, 12, 5688–5697 (2016)
      use param
      implicit none

      character(len = *), intent(in) :: filename
      integer  :: idof, iw, nmax
      real(dp) :: dw, w, tmp   ! angular freq. step
      real(dp) :: cutoff, scut ! for filtering freq. > fcut
      real(dp), dimension(:), allocatable :: fP
      complex(dp), dimension(:,:), allocatable :: Rf
      integer*8 :: plan ! for FFT

      print*, ''
      print*, '** Generating random forces **'

      ! number of points
      nmax = neq + nstep + 1

      ! ensure that nmax is even
      if (mod(nmax,2) /= 0) nmax = nmax + 1

      ! angular freq. step
      dw = TWOPI / (nmax * dt)
      scut = 3.d0 * dw

      allocate(fP(int(nmax / 2 - 1)))
      allocate(Rf(0:nmax-1,nrep))

      ! compute f_P function
      if (piqtbMode == 0) then
         call computefP0(fP, dw)
      else if (piqtbMode == 1) then
         call computefP1(fP, dw)
      else
         stop ('Error : wrong value of piqtbMode (= 0 or 1)')
      endif
      fP = dfloat(nrep) * KB * T * fP ! fP is now in energy units

      call FFT(Rf(:,1), Rf(:,1), plan, 'init') ! creates plan for FFT

      do idof = 1, nrep
         do iw = 1, int(nmax / 2 - 1)
            w = iw * dw
            ! cutoff for w > omegacut
            cutoff = 1.d0 / (dexp((w - omegacut) / scut) + 1.0d0)
            Rf(iw,idof) = dsqrt(fP(iw) * cutoff) * &
                           (gaussianRand(1.d0, 0.d0) * (1.0d0, 0.0d0) + &
                            gaussianRand(1.d0, 0.d0) * (0.0d0, 1.0d0))
         enddo
         ! to ensure that R(t) is real
         do iw = int(nmax / 2 + 1), nmax - 1
            Rf(iw,idof) = dconjg(Rf(nmax-iw,idof))
         enddo
         Rf(0,idof)  = (0.0d0, 0.0d0)
         Rf(nmax/2,idof)  = (0.0d0, 0.0d0)
         ! FFT of R(w) -> R(t)
         call FFT(Rf(:,idof), Rf(:,idof), plan, 'transform')
      enddo

      call FFT(Rf(:,1), Rf(:,1), plan, 'terminate')

      ! write random forces in file
      open(randfFileUnit, file = trim(filename), form = 'unformatted')
      ! normalisation : 1/sqrt(nmax) factor comes from FFT and 1/sqrt(2) to
      ! ensure that the random numbers (rk) used during the generation of
      ! the forces have a unit PSD (|rk|^2 = 1)
      ! see Brieuc et al., JCTC, 12, 5688–5697 (2016) - eq. (38) and (39)
      tmp = 1.d0 / dsqrt(2.d0 * nmax)
      do iw = 0, nmax - 1
         do idof = 1, nrep
            Rf(iw,idof) = Rf(iw,idof) * tmp
            if (piqtbMode == 1 .and. idof == 1) then
               ! centroid is classical
               Rf(iw,idof) = gaussianRand(dsqrt(KB*T), 0.d0) * (1.0d0, 0.0d0)
            endif
            write(randfFileUnit) real(realpart(Rf(iw,idof)))
         enddo
      enddo
      close(randfFileUnit)

      ! deallocation
      deallocate(Rf)
      deallocate(fP)

      print*, '** Random forces have been generated **'

   end subroutine generQTBRandF1

   subroutine generQTBRandF2(filename)
      ! generates forces associated with the QTB thermostat (ie. quantum
      ! Langevin thermostat)
      ! *********************************************
      ! if nrep = 1 : standard QTB thermostat
      ! *********************************************
      ! The random force is a gaussian colored noise with
      ! a Power Spectral Density (PSD) given by:
      ! 2 * m * gamma * theta(w, T) and
      ! theta(w, T) = hbar * w * ( 0.5 + 1 / (exp( hbar * w / kT ) - 1 ) )
      ! so that:
      ! <R(t)> = 0
      ! <R(t)R(t+tau) = FT(2 * m * gamma * theta(w,T)) (Wiener-Khinchin theorem)
      ! see Dammak et al, PRL, 103, 190601 (2009) for more details
      ! **********************************************
      ! if nrep > 1 : QTB-PIMD thermostat
      ! **********************************************
      ! The random force is also a stationnary gaussian colored noise but with
      ! a modified Power Spectral Density (PSD) given by:
      ! 2 * m * gamma * kappa(w, T) and
      ! kappa(w, T) = P * kT * f_P(hw/2kT) and f_P solution of
      ! eq. (13) for f_P^(0) and solution of eq. (16) for f_P^(1)
      ! so that:
      ! <R(t)> = 0
      ! <R(t)R(t+tau) = FT(2 * m * gamma * kappa(w,T)) (Wiener-Khinchin)
      ! see Brieuc, Dammak and Hayoun, JCTC, 12, 1351-1359 (2016)
      ! **********************************************
      ! the random force is generated here using the method proposed by Dammak
      ! and coworkers with a small modification : the random force is only
      ! generated for freq. < fcut and thus with a time step h = 1/(2*fcut).
      ! So random forces are not updated every timestep anymore but every
      ! mQTB = h / dt timesteps. The same rand. force applies for h timesteps.
      ! To account for this modification the target PSD has to be corrected.
      ! see Barrat and Rodney, J. Stat. Phys., 144, 679–689 (2011)
      use param
      implicit none

      character(len = *), intent(in) :: filename
      integer  :: idof, iw, nmax
      real(dp) :: w, dw      ! angular freq. and angular freq. step
      real(dp) :: correct, h ! correction of the spectrum to take into account
                             ! the time step h = mQTB * dt
      real(dp) :: tmp
      real(dp), dimension(:), allocatable :: fP
      complex(dp), dimension(:,:), allocatable :: Rf
      integer*8 :: plan ! for FFT

      print*, ''
      print*, '** Generating random forces **'

      ! number of points
      nmax = int(dfloat(neq + nstep + 1) / dfloat(mQTB)) + 1

      ! ensure that nmax is even
      if (mod(nmax,2) /= 0) nmax = nmax + 1

      ! angular freq. step - dw = 2*pi/(nmax*h) = 2*pi/(Ndt) = 2*wcut/nmax
      dw = 2.d0 * omegacut / nmax

      ! h = random force R(t) time step : h = mQTB * dt
      h = mQTB * dt

      ! allocation
      allocate(Rf(0:nmax-1,nrep))
      allocate(fP(int(nmax/2 - 1)))

      ! compute f_P function
      if (piqtbMode == 0) then
         call computefP0(fP, dw)
      else if (piqtbMode == 1) then
         call computefP1(fP, dw)
      else
         stop ('Error : wrong value of piqtbMode (= 0 or 1)')
      endif
      fP = dfloat(nrep) * KB * T * fP !fP is now in energy units

      call FFT(Rf(:,1), Rf(:,1), plan, 'init')

      do idof = 1, nrep
         do iw = 1, int(nmax / 2 - 1)
            w = iw * dw
            tmp = 0.5d0 * w * h
            correct = dsin(tmp) / tmp
            Rf(iw,idof) = dsqrt(fP(iw)) / correct * &
                        (gaussianRand(1.d0, 0.d0) * (1.0d0, 0.0d0) + &
                         gaussianRand(1.d0, 0.d0) * (0.0d0, 1.0d0))
         enddo
         ! to ensure that R(t) is real
         do iw = int(nmax / 2 + 1), nmax - 1
            Rf(iw,idof) = dconjg(Rf(nmax-iw,idof))
         enddo
         Rf(0,idof)  = (0.0d0, 0.0d0)
         Rf(nmax/2,idof)  = (0.0d0, 0.0d0)
         !FFT of R(w) -> R(t)
         call FFT(Rf(:,idof), Rf(:,idof), plan, 'transform')
      enddo

      call FFT(Rf(:,1), Rf(:,1), plan, 'terminate')

      ! write random forces in file
      open(randfFileUnit, file = trim(filename), form = 'unformatted')
      ! normalisation : 1/sqrt(nmax) factor comes from FFT and 1/sqrt(2)
      ! to ensure that the random numbers (rk) used during the generation of
      ! the forces have a unit PSD (|rk|^2 = 1)
      ! see Brieuc et al., JCTC, 12, 5688–5697 (2016) - eq. (38) and (39)
      tmp = 1.d0 / dsqrt(2.d0 * nmax)
      do iw = 0, nmax - 1
         do idof = 1, nrep
            Rf(iw,idof) = Rf(iw,idof) * tmp
            if (piqtbMode == 1 .and. idof == 1) then
               ! centroid is classical
               Rf(iw,idof) = gaussianRand(dsqrt(KB*T), 0.d0) * (1.0d0, 0.0d0)
            endif
            write(randfFileUnit) real(realpart(Rf(iw,idof)))
         enddo
      enddo
      close(randfFileUnit)

      ! deallocation
      deallocate(Rf)
      deallocate(fP)

      print*, '** Random forces have been generated **'

   end subroutine generQTBRandF2

   subroutine generQTBRandF3()
      ! generates forces associated with the QTB thermostat (ie. quantum
      ! Langevin thermostat)
      ! *********************************************
      ! if nrep = 1 : standard QTB thermostat
      ! *********************************************
      ! The random force is a gaussian colored noise with
      ! a Power Spectral Density (PSD) given by:
      ! 2 * m * gamma * theta(w, T) and
      ! theta(w, T) = hbar * w * ( 0.5 + 1 / (exp( hbar * w / kT ) - 1 ) )
      ! so that:
      ! <R(t)> = 0
      ! <R(t)R(t+tau) = FT(2 * m * gamma * theta(w,T)) (Wiener-Khinchin theorem)
      ! see Dammak et al, PRL, 103, 190601 (2009) for more details
      ! **********************************************
      ! if nrep > 1 : QTB-PIMD thermostat
      ! **********************************************
      ! The random force is also a stationnary gaussian colored noise but with
      ! a modified Power Spectral Density (PSD) given by:
      ! 2 * m * gamma * kappa(w, T) and
      ! kappa(w, T) = P * kT * f_P(hw/2kT) and f_P solution of
      ! eq. (13) for f_P^(0) and solution of eq. (16) for f_P^(1)
      ! so that:
      ! <R(t)> = 0
      ! <R(t)R(t+tau) = FT(2 * m * gamma * kappa(w,T)) (Wiener-Khinchin)
      ! see Brieuc, Dammak and Hayoun, JCTC, 12, 1351-1359 (2016)
      ! **********************************************
      ! the random force is generated "on-the-fly" using the method proposed by
      ! Barrat and Rodney. The random force is generated for freq. < fcut
      ! and thus with a time step h = 1/(2*fcut). So random forces are not
      ! updated every timestep but every mQTB = h / dt timesteps.
      ! The same rand. force applies for h timesteps. To account for this
      ! modification the target PSD has to be corrected. This subroutine
      ! initialize the values of the filter H and of the first random numbers r.
      ! The forces are then generated at every mQTB steps by computing
      ! the convolution product.
      ! see Barrat and Rodney, J. Stat. Phys., 144, 679–689 (2011)
      use param
      implicit none

      integer  :: irep, iw, iiw, it, nmax
      real(dp) :: w, dw       ! angular freq. and angular freq. step
      real(dp) :: correct, h  ! correction of the spectrum to take into account
                              ! the time step h = mQTB * dt
      real(dp) :: tmp, tt
      real(dp), dimension(:), allocatable :: Hw ! H(w)
      real(dp), dimension(:), allocatable :: fP

      ! ensure that nmax is even
      ! if (mod(nQTB,2) /= 0) nQTB = nQTB + 1

      print*, ''
      print*, '** Initializing random forces **'

      ! number of points
      nmax = nQTB

      ! h = random force R(t) time step : h = mQTB * dt
      h = mQTB * dt

      ! angular freq. step
      dw = TWOPI / (nmax * h)

      allocate(fP(nmax/2))
      allocate(Hw(0:nmax - 1))

      ! compute f_P function
      if (piqtbMode == 0) then
         call computefP0(fP, dw)
      else if (piqtbMode == 1) then
         call computefP1(fP, dw)
      else
         stop ('Error : wrong value of piqtbMode (= 0 or 1)')
      endif
      fP = dfloat(nrep) * KB * T * fP !fP is now in energy units

      ! compute the "filter" in frequency space H(w)
      do iw = 0, nmax/2 - 1
         iiw = iw - nmax / 2
         w = iiw * dw
         tmp = 0.5d0 * w * h
         correct = dsin(tmp) / tmp
         Hw(iw) = dsqrt(fP(-iiw)) / correct
         ! write(200,*) iw, au2THz(iw*dw/TWOPI), Hw(iw)
      enddo
      Hw(nmax/2) = 0.d0
      ! write(200,*) nmax/2, au2THz(nmax/2*dw/TWOPI), Hw(nmax/2)
      do iw = nmax/2 + 1, nmax - 1
         iiw = iw - nmax / 2
         w = iiw * dw
         tmp = 0.5d0 * w * h
         correct = dsin(tmp) / tmp
         Hw(iw) = dsqrt(fP(iiw)) / correct
         ! write(200,*) iw, au2THz(iw*dw/TWOPI), Hw(iw)
      enddo

      ! compute the "filter" in time space H(t)
      do it = 0, nmax - 1
         HQTB(it) = 0.d0
         tt = (it - nmax / 2)
         do iw = 0, nmax - 1
            w = (iw - nmax / 2) * TWOPI / nmax
            HQTB(it) = HQTB(it) + Hw(iw) * dcos(w*tt)
         enddo
         HQTB(it) = HQTB(it) / nmax
         ! initialize random numbers
         do irep = 1, nrep
            rQTB(irep,it) = gaussianRand(1.d0,0.d0)
         enddo
      enddo

      deallocate(Hw)
      deallocate(fP)

      print*, '** Random forces are initialized **'

   end subroutine generQTBRandF3

   subroutine computefP0(fP, dw)
      ! compute the solution F_P^(0) of equation (13) in
      ! Brieuc, Dammak and Hayoun, JCTC, 12, 1351-1359 (2016)
      use param
      implicit none

      real(dp), dimension(:), intent(out) :: fP ! the solution
      real(dp), intent(in) :: dw    ! angular freq. step
      real(dp) :: xmin, xmax, dx    ! min value and step for x = sqrt(P)*u (F_P)
      real(dp) :: malpha            ! mixing parameter alpha
      integer  :: niter             ! nb of iterations
      real(dp), dimension(:), allocatable   :: x, x2, h, fP1 ! x = sqrt(P)*u
      real(dp), dimension(:,:), allocatable :: xk, xk2, fPxk ! "xk = uk(x)"
      integer, dimension(:,:), allocatable  :: kk       ! fPxk = F_P(xk*sqrt(P))
      integer  :: i, j, k, n, nx
      real(dp) :: w, HBokT, tmp, OSQNREP, fPrev, err
      real(dp) :: aa,bb,r2 ! for linear regression
      real(dp) :: x1, dx1

      n = size(fP)
      HBokT = HB / (KB * T)

      ! P = 1 : standard QTB
      ! fP = theta(w, T) / kT
      if (nrep == 1) then
         do k = 1, n
            w = dfloat(k) * dw
            tmp = HBokT * w
            fP(k) = tmp * (0.5d0 + 1.d0 / (exp(tmp) - 1.d0))
         enddo

      ! P > 2 : QTB-PIMD
      else if (nrep > 1) then
         print*, 'computing f_P (0) function'
         ! **** initialization ****
         dx1 = dsqrt(dfloat(nrep)) * 0.5d0 * HBokT * dw ! x = u*sqrt(P)
         xmin = 1.0d-7
         dx = 0.05d0
         xmax = 10000.d0
         nx = int((xmax - xmin)/dx) + 1
         nx = nx + nx/5  !add 20% points to avoid problem at the end
                         !of the interval
         OSQNREP = 1.d0 / dsqrt(dfloat(nrep))
         malpha = ONREP ! alpha = 1/P
         niter = 30     ! 30 iterations are enough to converge

         print*, 'dw =', dw
         print*, 'dx =', dx
         print*, 'xmin:', xmin
         print*, 'xmax =', xmax
         print*, 'nx, n =', nx, n

         ! allocation
         allocate(x(nx))
         allocate(x2(nx))
         allocate(h(nx))
         allocate(fP1(nx))
         allocate(xk(1:nrep-1,nx))
         allocate(xk2(1:nrep-1,nx))
         allocate(kk(1:nrep-1,nx))
         allocate(fPxk(1:nrep-1,nx))

         ! initialize F_P(x) = f_P(x/sqrt(P))
         ! fP -> F_P(x) = 1/P * h(x/P)
         ! fPxk -> F_P(xk*sqrt(P)) = 1/P * h(xk/sqrt(P))
         do j = 1, nx
            x(j) = xmin + (j - 1) * dx   ! x
            x2(j) = x(j)**2              ! x^2
            h(j) = x(j)/dtanh(x(j))      ! h(x)
            if (x(j) <= 1.d-10) h(j)=1.d0 ! to avoid division by 0
            fP1(j) = ONREP*x(j)*ONREP/dtanh(x(j)*ONREP) ! F_P(x)
            if (x(j)*ONREP <= 1.d-10) fP1(j)=ONREP ! to avoid division by 0
            do k = 1, nrep-1
               xk2(k,j) = ONREP*x2(j)+dfloat(nrep)*dsin(dfloat(k)*PI*ONREP)**2
               xk(k,j) = dsqrt(nrep*xk2(k,j)) ! sqrt(P)*xk
               kk(k,j) = nint((xk(k,j) - xmin) / dx) + 1 ! k <-> sqrt(P)*xk
               fpxk(k,j) = ONREP*xk(k,j)*ONREP/dtanh(xk(k,j)*ONREP)! F_P(sqrt(P)*xk)
               if(xk(k,j)*ONREP <= 1.d-10) fPxk(k,j) = ONREP
            enddo
         enddo

         ! **** resolution ****
         ! compute F_P(x)
         do i = 1, niter
            err = 0.d0
            do j = 1, nx
               tmp = 0.d0
               do k = 1, nrep-1
                  tmp = tmp + fpxk(k,j) * x2(j) / xk2(k,j)
               enddo
               fPrev = fP1(j)
               fP1(j) = malpha * ONREP * (h(j) - tmp) + (1.d0 - malpha) * fP1(j)
               if (j <= n) err = err + dabs(1-fP1(j)/fPrev) ! compute "errors"
            enddo
            err = err / dfloat(n)

            ! Linear regression on the last 20% of the F_P function
            call linReg(fP1(8*nx/10:nx),x(8*nx/10:nx),aa,bb,r2)

            ! compute the new F_P(xk*sqrt(P))
            ! through linear interpolation
            ! or linear extrapolation if outside of the range
            do j = 1, nx
               do k = 1, nrep-1
                  if (kk(k,j) < nx) then
                     ! linear interpolation
                     fPxk(k,j) = fP1(kk(k,j))+(fP1(kk(k,j)+1)-fP1(kk(k,j)))/dx &
                     &           * (xk(k,j) - x(kk(k,j)))
                  else
                     ! linear extrapolation
                     fPxk(k,j) = aa * xk(k,j) + bb
                  endif
               enddo
            enddo
         enddo

         ! **** tests ****
         open(100,file='fP_function.res')
         ! "average error"
         if (err > 1.d-4) then
            print*, 'probable error during the generation of F_P function :'
            print*, 'average "error"', err
         endif
         ! value at zero freq.
         print*, 'F_P at zero freq. - theoretical (= 1/P):', ONREP
         print*, 'F_P at zero freq. - calculated:', fP1(1)

         ! slope at high freq.
         if (dabs(1.d0-aa/(ONREP*ONREP)) > 1e-2) then
            print*, 'probable error during the generation of F_P function :'
            print*, dabs(1.d0-aa/(ONREP*ONREP))
            print*, 'slope of F_P at high freq. - theoretical:', ONREP*OSQNREP
            print*, 'slope of F_P at high freq. - calculated:', aa/OSQNREP
         else
            print*, 'slope of F_P at high freq. - theoretical:', ONREP*OSQNREP
            print*, 'slope of F_P at high freq. - calculated:', aa/OSQNREP
         endif

         write(100,*) '# average error during computation', err
         write(100,*) '# slope of F_P at high freq. - theoretical:', ONREP*OSQNREP
         write(100,*) '# slope of F_P at high freq. - calculated:', aa/OSQNREP
         write(100,* )'# F_P at zero freq. - theoretical (= 1/P):', ONREP
         write(100,* )'# F_P at zero freq. - calculated:', fP1(1)

         write(100,* )''

         ! **** write solution ****
         ! Linear regression on the last 20% of the F_P function
         call linReg(fP1(8*nx/10:nx/10),x(8*nx/10:nx/10),aa,bb,r2)
         do j = 1, n
            x1 = j * dx1
            k = nint((x1 - xmin) / dx) + 1
            if (k > nx) then
               fP(j) = aa * x1 + bb
            else if (k <= 0) then
               stop ('error in fP computation x < xmin')
            else
               if (x1 > x(k)) then
                  fP(j) = fP1(k) + (fP1(k+1)-fP1(k))/dx * (x1-x(k))
               else
                  fP(j) = fP1(k) + (fP1(k)-fP1(k-1))/dx * (x1-x(k))
               endif
            endif
         enddo

         do j = 1, nx
            write(100,*) j, au2THz(j*dw/TWOPI), x(j)*OSQNREP, fP1(j)
         enddo

         ! deallocation
         deallocate(x)
         deallocate(x2)
         deallocate(h)
         deallocate(xk)
         deallocate(xk2)
         deallocate(kk)
         deallocate(fPxk)
         deallocate(fP1)

         print*, 'fP function has been computed'
         close(100)

      endif

   end subroutine computefP0

   subroutine computefP1(fP, dw)
      ! compute the solution F_P^(1) of equation (16) in
      ! Brieuc, Dammak and Hayoun, JCTC, 12, 1351-1359 (2016)
      use param
      implicit none

      real(dp), dimension(:), intent(out) :: fP ! the solution
      real(dp), intent(in) :: dw    ! angular freq. step
      real(dp) :: xmin, xmax, dx    ! min value and step for x = sqrt(P)*u (F_P)
      real(dp) :: malpha            ! mixing parameter alpha
      integer  :: niter             ! nb of iterations
      real(dp), dimension(:), allocatable   :: x, x2, h, fP1 ! x = sqrt(P)*u
      real(dp), dimension(:,:), allocatable :: xk, xk2, fPxk ! "xk = uk(x)"
      integer, dimension(:,:), allocatable  :: kk       ! fPxk = F_P(xk*sqrt(P))
      integer  :: i, j, k, n, nx
      real(dp) :: w, HBokT, tmp, tmp1, OSQNREP, ONREP1, fPrev, err
      real(dp) :: aa,bb,r2 ! for linear regression
      real(dp) :: dx1

      n = size(fP)
      HBokT = HB / (KB * T)

      ! P = 1 : standard QTB
      ! fP = theta(w, T) / kT
      if (nrep == 1) then
         do k = 1, n
            w = dfloat(k) * dw
            tmp = HBokT * w
            fP(k) = tmp * (0.5d0 + 1.d0 / (exp(tmp) - 1.d0))
         enddo

      ! P > 2 : QTB-PIMD
      else if (nrep > 1) then
         print*, 'computing f_P (1) function'
         ! **** initialization ****
         dx1 = dsqrt(dfloat(nrep)) * 0.5d0 * HBokT * dw ! x = u*sqrt(P)
         xmin = 1.0d-7
         dx = 0.05d0
         xmax = 10000.d0
         nx = int((xmax - xmin)/dx) + 1
         nx = nx + nx/5  !add 20% points to avoid problem at the end
                         !of the interval
         OSQNREP = 1.d0 / dsqrt(dfloat(nrep))
         ONREP1 = 1.d0 / dfloat(nrep - 1)
         malpha = ONREP ! alpha = 1/P
         niter = 30     ! 30 iterations are enough to converge

         print*, 'dw=', dw
         print*, 'dx=', dx
         print*, 'xmax=', xmax
         print*, 'nx, n, mu=', nx, n

         ! allocation
         allocate(x(nx))
         allocate(x2(nx))
         allocate(h(nx))
         allocate(fP1(nx))
         allocate(xk(2:nrep-1,nx))
         allocate(xk2(2:nrep-1,nx))
         allocate(kk(2:nrep-1,nx))
         allocate(fPxk(2:nrep-1,nx))

         ! initialize F_P(x) = f_P(x/sqrt(P))
         ! fP -> F_P(x) = [h(x/P)-1/P]/(P-1)
         ! fPxk -> F_P(sqrt(P*xk**2-(P*sin(pi/P))**2)
         do j = 1, nx
            x(j) = xmin + (j - 1) * dx ! x
            x2(j) = x(j)**2            ! x^2
            h(j) = x(j)/dtanh(x(j))    ! h(x)
            if(x(j) <= 1.d-10) h(j)=1.d0
            fP1(j) = ONREP1*ONREP * (x(j)/dtanh(x(j)*ONREP) - 1.d0) ! F_P(x)
            if(x(j)*ONREP <= 1.d-10) fp1(j)=ONREP
            do k = 2, nrep-1
               xk2(k,j) = ONREP*x2(j)+dfloat(nrep)*dsin(dfloat(k)*PI*ONREP)**2
               xk(k,j) = dsqrt(nrep*xk2(k,j) - (dfloat(nrep)*dsin(PI*ONREP))**2)
               kk(k,j) = nint((xk(k,j) - xmin) / dx) + 1
               fpxk(k,j) = ONREP1*ONREP * (xk(k,j)/dtanh(xk(k,j)*ONREP) - 1.d0)
               if(xk(k,j)*ONREP <= 1.d-10) fPxk(k,j) = ONREP
            enddo
         enddo

         ! **** resolution ****
         ! compute F_P(x)
         do i = 1, niter
            err = 0.d0
            do j = 1, nx
               tmp = 0.d0
               do k = 2, nrep-1
                  tmp = tmp + fpxk(k,j) * x2(j) / xk2(k,j)
               enddo
               fPrev = fP1(j)
               tmp1 = ONREP + dfloat(nrep) * (dsin(PI*ONREP) / x(j))**2
               fP1(j) = malpha*tmp1 * (h(j)-1.d0-tmp) + (1.d0 - malpha) * fP1(j)
               if (j <= n) err = err + dabs(1 - fP1(j)/fPrev) ! compute "errors"
            enddo
            err = err / dfloat(n)

            ! Linear regression on the last 10% of the F_P function
            call linReg(fP1(8*nx/10:9*nx/10),x(8*nx/10:9*nx/10),aa,bb,r2)

            ! compute the new F_P(sqrt(P*xk**2-(P*sin(pi/P))**2)
            ! through linear interpolation
            ! or linear extrapolation if outside of the range
            do j = 1, nx
               do k = 2, nrep-1
                  if (kk(k,j) < nx) then
                     ! linear interpolation
                     fPxk(k,j) = fP1(kk(k,j))+(fP1(kk(k,j)+1)-fP1(kk(k,j)))/dx &
                     &           * (xk(k,j) - x(kk(k,j)))
                  else
                     ! linear extrapolation
                     fPxk(k,j) = aa * xk(k,j) + bb
                  endif
               enddo
            enddo
         enddo

         ! **** tests ****
         open(100,file='fP_function.res')
         ! "average error"
         if (err > 1.d-4) then
            print*, 'probable error during the generation of F_P function :'
            print*, 'average "error"', err
         endif
         ! value at zero freq.
         ! print*, 'F_P at zero freq. - theoretical (= 1/P):', ONREP
         ! print*, 'F_P at zero freq. - calculated:', fP1(1)

         ! slope at high freq.
         if (dabs(1.d0-aa/(ONREP1*ONREP)) > 1e-2) then
            print*, 'probable error during the generation of F_P function :'
            print*, dabs(1.d0-aa/(ONREP1*ONREP))
            print*, 'slope of F_P at high freq. - theoretical:', ONREP1*OSQNREP
            print*, 'slope of F_P at high freq. - calculated:', aa/OSQNREP
         else
            print*, 'slope of F_P at high freq. - theoretical:', ONREP1*OSQNREP
            print*, 'slope of F_P at high freq. - calculated:', aa/OSQNREP
         endif

         write(100,*) '# average error during computation', err
         write(100,*) '# slope of F_P at high freq. - theoretical:', ONREP1*OSQNREP
         write(100,*) '# slope of F_P at high freq. - calculated:', aa/OSQNREP

         ! **** write solution ****
         ! Linear regression on the last 20% of the F_P function
         call linReg(fP1(8*nx/10:nx/10),x(8*nx/10:nx/10),aa,bb,r2)
         do j = 1, n
            tmp = dfloat(j) * dx1 * OSQNREP
            tmp = dfloat(nrep) * (tmp**2 - dfloat(nrep) * dsin(PI*ONREP)**2)
            if (tmp < 0.d0) then
               fP(j) = 0.0d0
            else
               tmp = dsqrt(tmp)
               k = nint((tmp - xmin) / dx)
               if(k <= 0) then
                  fP(j) = fP1(1)
               else if (k > nx) then
                  fP(j) = aa * tmp + bb
               else
                  if (tmp > x(k)) then
                     fP(j) = fP1(k) + (fP1(k+1)-fP1(k))/dx * (tmp - x(k))
                  else
                     fP(j) = fP1(k) + (fP1(k)-fP1(k-1))/dx * (tmp - x(k))
                  endif
               endif
            endif
         enddo

         do j = 1, nx
            write(100,*) j, au2THz(j*dw/TWOPI), x(j)*OSQNREP, fP1(j)
         enddo

         ! deallocation
         deallocate(x)
         deallocate(x2)
         deallocate(h)
         deallocate(xk)
         deallocate(xk2)
         deallocate(kk)
         deallocate(fPxk)
         deallocate(fP1)

         print*, 'fP function has been computed'
         close(100)

      endif

   end subroutine computefP1

   subroutine FFT(x_in, x_out, plan, mode)
      implicit none

      complex*16, dimension(:), intent(in)  :: x_in
      complex*16, dimension(:), intent(out) :: x_out
      character(len=*), intent(in) :: mode
      integer*8, intent(inout) :: plan ! required for FFTW

      integer :: n, ntest

      n = size(x_in)
      ntest = size(x_out)
      if (n /= ntest) stop ('Error with FFT during QTB random force generation')

      if (trim(mode) == 'init') then
         call dfftw_plan_dft_1d(plan,n,x_in,x_out,'FFTW_FORWARD','FTTW_ESTIMATE')
      else if (trim(mode) == 'transform') then
         ! Forward FFT using FFTW3 library
         call dfftw_execute_dft(plan,x_in,x_out)
      else if (trim(mode) == 'terminate') then
         call dfftw_destroy_plan(plan)
      endif

   end subroutine FFT

   real(dp) function gaussianRand(stdDev, av)
      ! generates random numbers in a gaussian distribution with
      ! average 'av' and standard deviation 'std'
      ! Based on the Box-Muller transformation
      use param
      implicit none

      real(dp), intent(in) :: stdDev
      real(dp), intent(in) :: av
      real(dp) :: ZBQLU01 ! generates uniform random number in [0,1]
      real(dp) :: random

      random = dcos(TWOPI * ZBQLU01()) * dsqrt(-2.d0 * dlog(ZBQLU01()))
      gaussianRand = stdDev * random + av

   end function gaussianRand

   subroutine linReg(y,x,a,b,r2)
      ! simple linear regression of y
      ! so that y(x) = a*x + b
      implicit none

      real(dp), dimension(:), intent(in) :: y, x
      real(dp), intent(out) :: a, b
      real(dp), intent(out) :: r2 ! coefficient of determination r^2
      integer :: i, n
      real(dp) :: xav, yav   ! x and y average value
      real(dp) :: xycov      ! xy covariance
      real(dp) :: xvar, yvar ! x and y variance

      n = size(y)
      if(size(x) /= n) then
         stop ('error in linReg subroutine: x and y don''t have the same size')
      endif

      xav = 0.d0
      yav = 0.d0
      xycov = 0.d0
      xvar = 0.d0
      yvar = 0.d0

      do i = 1, n
         xav = xav + x(i)
         yav = yav + y(i)
         xycov = xycov + x(i)*y(i)
         xvar = xvar + x(i)**2
         yvar = yvar + y(i)**2
      enddo

      xav = xav / dfloat(n)
      yav = yav / dfloat(n)
      xycov = xycov / dfloat(n)
      xycov = xycov - xav*yav
      xvar = xvar / dfloat(n)
      xvar = xvar - xav**2
      yvar = yvar / dfloat(n)
      yvar = yvar - yav**2

      a = xycov / xvar
      b = yav - a * xav

      r2 = xycov / dsqrt(xvar * yvar)

      if (r2 < 0.9d0) write(6,'(a85, e10.3)') 'error in the linear regression &
      & part of the f_P function generation is probable : r^2 =', r2

   end subroutine linReg

end module thermostats
