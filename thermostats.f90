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

      character(len=*), intent(in) :: mode

      if (mode == 'initialize') then
         select case(therm)
            case('langevin')
               std = dsqrt(2.D0 * gam * M * KB * T / dt)
            case('qtb')
               ! Random force generation : Standard Dammak et al
               if (noiseGenMode == 1) then
                  ! new random force at every steps
                  mQTB = 1
                  std = dsqrt(2.D0 * gam * M / (mQTb * dt))
                  call generQTBRandF1('qtbRandF.tmp')
                  open(randfFileUnit, file ='qtbRandF.tmp', form='unformatted')
               ! Random force generation : Modified Dammak et al
               else if (noiseGenMode == 2) then
                  ! new random force at every mQTB = 1/(2*fmax*dt) steps
                  mQTB = nint(1.d0 / (2.d0 * fmax * dt))
                  if (mQTB == 0) mQTB = 1
                  print*, 'mQTB =', mQTB
                  std = dsqrt(2.D0 * gam * M / (mQTb * dt))
                  call generQTBRandF2('qtbRandF.tmp')
                  open(randfFileUnit, file ='qtbRandF.tmp', form='unformatted')
               ! Random force generation : Barrat and Rodney
               else if (noiseGenMode == 3) then
                  ! new random force at every mQTB = 1/(2*fmax*dt) steps
                  mQTB = nint(1.d0 / (2.d0 * fmax * dt))
                  if (mQTB == 0) mQTB = 1
                  print*, 'mQTB =', mQTB
                  std = dsqrt(2.D0 * gam * M / (mQTb * dt))
                  ! ensure that nQTB is even
                  if (mod(nQTB,2) /= 0) nQTB = nQTB + 1
                  allocate(HQTB(0:nQTB-1)) !"filter" (PSD) in the time domain
                  allocate(rQTB(nrep,0:nQTB-1)) !random numbers
                  call generQTBRandF3()

               endif
         end select
      else if (mode == 'terminate') then
         select case(therm)
            case('qtb')
               if (noiseGenMode == 1 .or. noiseGenMode == 2) then
                  close(randfFileUnit)
               else if (noiseGenMode == 3) then
                  deallocate(HQTB)
                  deallocate(rQTB)
               endif
               call system('rm -f qtbRandF.tmp')
         end select
      endif

   end subroutine thermostating

   real(dp) function langevin(v)
      ! compute the forces associated with the standard (ie. classical)
      ! Langevin thermostat
      ! ******************************************
      ! The random fore is a gaussian white noise:
      ! <R(t)> = 0
      ! <R(t)R(t+tau)> = 2 * m * gam * kT * delta(tau)
      use param
      implicit none

      real(dp), intent(in) :: v

      R = gaussianRand(std, 0.d0)
      langevin = R - m * gam * v

   end function langevin

   real(dp) function qtb(v, irep, n)
      ! compute forces associated with the QTB thermostat (ie. quantum
      ! Langevin thermostat)
      use param
      implicit none

      real(dp), intent(in) :: v ! velocity
      integer, intent(in)  :: n ! step number
      integer, intent(in)  :: irep ! replica nb
      integer :: i
      real :: R_sngl

      !new QTB random force every mQTB step
      if (mod(n, mQTB) == 0) then
         if (noiseGenMode == 1 .or. noiseGenMode == 2) then
            read(randfFileUnit) R_sngl
            R = R_sngl
         else if (noiseGenMode == 3) then
            ! compute random force through the convolution product
            R = 0.d0
            do i = 0, nQTB-2
               R = R + HQTB(i) * rQTB(irep, i)
               rQTB(irep, i) = rQTB(irep, i + 1)
            enddo
            R = R + HQTB(nQTB - 1) * rQTB(irep, nQTB - 1)
            rQTB(irep, nQTB - 1) = gaussianRand(1.d0,0.d0)
         endif
      endif

      qtb = R - m * gam * v

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
      ! theta(w, T) = 0.5 * hbar * w * ( 1 + 1 / (exp( hbar * w / kT ) - 1 ) )
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
      real(dp) :: dw !angular freq. step
      real(dp) :: cutoff, scut ! for filtering freq. > fmax
      real(dp), dimension(:), allocatable :: fP
      complex(dp), dimension(:,:), allocatable :: Rw, Rt

      print*, ' ** Generating random forces **'

      ! number of points
      nmax = neq + nstep + 1

      ! ensure that nmax is even
      if (mod(nmax,2) /= 0) nmax = nmax + 1

      ! angular freg. step
      dw = TWOPI / (nmax * dt)
      scut = 3.d0 * dw

      ! allocation
      allocate(Rw(nrep,0:nmax-1))
      allocate(Rt(nrep,0:nmax-1))
      allocate(fP(int(nmax / 2 - 1)))

      ! compute f_P function
      if (piqtbMode == 0) then
         call computefP0(fP, dw)
      else if (piqtbMode == 1) then
         call computefP1(fP, dw)
      else
         stop('Error : wrong value of piqtbMode (= 0 or 1)')
      endif
      fP = dfloat(nrep) * KB * T * fP !fP is now in au energy units (Ha)

      do idof = 1, nrep
         do iw = 1, int(nmax / 2 - 1)
            ! cutoff for w > omegamax
            cutoff = 1.d0 / (dexp((iw*dw - omegamax) / scut) + 1.0d0)
            Rw(idof, iw) = dsqrt(fP(iw) * cutoff) * &
            &(gaussianRand(1.d0, 0.d0) * (1.0d0, 0.0d0) + &
            & gaussianRand(1.d0, 0.d0) * (0.0d0, 1.0d0))
         enddo
         ! to ensure that R(t) is real
         do iw = int(nmax / 2 + 1), nmax - 1
            Rw(idof, iw) = dconjg(Rw(idof, nmax - iw))
         enddo
         Rw(idof, 0)  = (0.0d0, 0.0d0)
         Rw(idof, nmax / 2)  = (0.0d0, 0.0d0)
         !FFT of R(w) -> R(t)
         call FFT(Rw(idof,:), Rt(idof,:))
      enddo

      ! write random forces in file
      open(randfFileUnit, file = trim(filename), form = 'unformatted')
      do iw = 0, nmax - 1
         do idof = 1, nrep
            Rt(idof,iw) = Rt(idof,iw) * std / dsqrt(2.d0 * nmax)
            write(randfFileUnit) real(realpart(Rt(idof,iw)))
         enddo
      enddo
      close(randfFileUnit)

      ! deallocation
      deallocate(Rw)
      deallocate(Rt)
      deallocate(fP)

      print*, ' ** Random forces have been generated **'

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
      ! theta(w, T) = 0.5 * hbar * w * ( 1 + 1 / (exp( hbar * w / kT ) - 1 ) )
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
      ! generated for freq. < fmax and thus with a time step h = 1/(2*fmax).
      ! So random forces are not updated every timestep anymore but every
      ! mQTB = h / dt timesteps. The same rand. force applies for h timesteps.
      ! To account for this modification the target PSD as to be corrected.
      ! see Barrat and Rodney, J. Stat. Phys., 144, 679–689 (2011)
      use param
      implicit none

      character(len = *), intent(in) :: filename
      integer  :: idof, iw, nmax
      real(dp) :: w, dw !angular freq. and angular freq. step
      real(dp) :: correct, h ! correction of the spectrum to take into account
                             ! the time step h = mQTB * dt
      real(dp) :: tmp
      real(dp), dimension(:), allocatable :: fP
      complex(dp), dimension(:,:), allocatable :: Rw, Rt

      print*, ' ** Generating random forces **'

      ! number of points
      nmax = int(dfloat(neq + nstep + 1) / dfloat(mQTB)) + 1

      ! ensure that nmax is even
      if (mod(nmax,2) /= 0) nmax = nmax + 1

      ! angular freg. step
      dw = 2.d0 * omegamax / nmax

      ! h = ranom force R(t) time step : h = mQTB * dt
      h = mQTB * dt

      ! allocation
      allocate(Rw(nrep,0:nmax-1))
      allocate(Rt(nrep,0:nmax-1))
      allocate(fP(int(nmax / 2 - 1)))

      ! compute f_P function
      if (piqtbMode == 0) then
         call computefP0(fP, dw)
      else if (piqtbMode == 1) then
         call computefP1(fP, dw)
      else
         stop('Error : wrong value of piqtbMode (= 0 or 1)')
      endif
      fP = dfloat(nrep) * KB * T * fP !fP is now in au energy units (Ha)

      do idof = 1, nrep
         do iw = 1, int(nmax / 2 - 1)
            w = iw * dw
            tmp = 0.5d0 * w * h
            correct = dsin(tmp) / tmp
            Rw(idof, iw) = dsqrt(fP(iw)) / correct * &
            &(gaussianRand(1.d0, 0.d0) * (1.0d0, 0.0d0) + &
            & gaussianRand(1.d0, 0.d0) * (0.0d0, 1.0d0))
         enddo
         ! to ensure that R(t) is real
         do iw = int(nmax / 2 + 1), nmax - 1
            Rw(idof, iw) = dconjg(Rw(idof, nmax - iw))
         enddo
         Rw(idof, 0)  = (0.0d0, 0.0d0)
         Rw(idof, nmax / 2)  = (0.0d0, 0.0d0)
         !FFT of R(w) -> R(t)
         call FFT(Rw(idof,:), Rt(idof,:))
      enddo

      ! write random forces in file
      open(randfFileUnit, file = trim(filename), form = 'unformatted')
      do iw = 0, nmax - 1
         do idof = 1, nrep
            Rt(idof,iw) = Rt(idof,iw) * std / dsqrt(2.d0 * nmax)
            write(randfFileUnit) real(realpart(Rt(idof,iw)))
         enddo
      enddo
      close(randfFileUnit)

      ! deallocation
      deallocate(Rw)
      deallocate(Rt)
      deallocate(fP)

      print*, ' ** Random forces have been generated **'

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
      ! theta(w, T) = 0.5 * hbar * w * ( 1 + 1 / (exp( hbar * w / kT ) - 1 ) )
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
      ! Barrat and Rodney. The random force is generated for freq. < fmax
      ! and thus with a time step h = 1/(2*fmax). So random forces are not
      ! updated every timestep but every mQTB = h / dt timesteps.
      ! The same rand. force applies for h timesteps. To account for this
      ! modification the target PSD as to be corrected. This subroutine
      ! initialize the values of the filter and of the first random numbers.
      ! The forces are then generated at every mQTB steps by computing
      ! the convolution product.
      ! see Barrat and Rodney, J. Stat. Phys., 144, 679–689 (2011)
      use param
      implicit none

      integer  :: idof, iw, nmax
      real(dp) :: w, dw !angular freq. and angular freq. step
      real(dp) :: correct, h ! correction of the spectrum to take into account
                             ! the time step h = mQTB * dt
      real(dp) :: tmp
      complex(dp), dimension(:), allocatable :: Hw, Ht
      real(dp), dimension(:), allocatable :: fP

      ! ensure that nmax is even
      ! if (mod(nQTB,2) /= 0) nQTB = nQTB + 1

      print*, ' ** Initializing random forces **'

      ! number of points
      nmax = nQTB

      ! angular freg. step
      dw = 2.d0 * omegamax / nmax

      ! h = ranom force R(t) time step : h = mQTB * dt
      h = mQTB * dt

      allocate(fP(int(nmax / 2 - 1)))
      allocate(Hw(0:nmax - 1))
      allocate(Ht(0:nmax - 1))

      ! compute f_P function
      if (piqtbMode == 0) then
         call computefP0(fP, dw)
      else if (piqtbMode == 1) then
         call computefP1(fP, dw)
      else
         stop('Error : wrong value of piqtbMode (= 0 or 1)')
      endif
      fP = dfloat(nrep) * KB * T * fP !fP is now in au energy units (Ha)

      ! compute the "filter" H(w)
      do iw = 1, int(nmax / 2 - 1)
         w = iw * dw
         tmp = 0.5d0 * w * h
         correct = dsin(tmp) / tmp
         Hw(iw) = dsqrt(fP(iw)) / correct * (1.d0, 0.d0)
      enddo
      ! to ensure that H(t) is real
      do iw = int(nmax / 2 + 1), nmax - 1
         Hw(iw) = dconjg(Hw(nmax - iw))
      enddo
      Hw(0)  = (0.0d0, 0.0d0)
      Hw(nmax / 2)  = (0.0d0, 0.0d0)
      ! FFT of H(w) -> H(t)
      call FFT(Hw, Ht)
      do iw = 0, nmax - 1
         HQTB(iw) = realpart(Ht(iw)) * std / dfloat(nmax)
         do idof = 1, nrep
            rQTB(idof, iw) = gaussianRand(1.d0,0.d0)
         enddo
      enddo

      deallocate(Hw)
      deallocate(Ht)

      print*, ' ** Random forces are initialized **'

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
         dx = dsqrt(dfloat(nrep)) * 0.5d0 * HBokT * dw ! x = u*sqrt(P)
         xmin = dx
         xmax = dsqrt((n*dx)**2 + dfloat(nrep)**2) ! sqrt(P)*x_kmax
         if (xmax < 50.d0) xmax = 50.d0    ! xmax >= 50 to avoid errors
         nx = int((xmax - xmin)/dx) + 1
         nx = nx + nx/5  !add 20% points to avoid problem at the end
                         !of the interval
         OSQNREP = 1.d0 / dsqrt(dfloat(nrep))
         malpha = ONREP ! alpha = 1/P
         niter = 30     ! 30 iterations are enough to converge

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
               xk(k,j) = dsqrt(nrep*xk2(k,j)) !sqrt(P)*xk
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

            ! Linear regression on the last 10% of the F_P function
            call linReg(fP1(9*n/10:n),x(9*n/10:n),aa,bb,r2)

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
         ! "average error"
         if (err > 1.d-6) then
            print*, 'error during the generation of the F_P function is highly &
            &        probable :'
            print*, 'average "error"', err
         endif
         ! value at zero freq.
         if (dabs(1.d0-fP1(1)/ONREP) > 1e-6) then
            print*, 'error during the generation of the F_P function is highly &
            &        probable :'
            print*, 'F_P at zero freq. - theoretical (= 1/P):', ONREP
            print*, 'F_P at zero freq. - calculated:', fP1(1)
         endif
         ! slope at high freq.
         if (dabs(1.d0-aa/(ONREP*ONREP)) > 5e-3) then
            print*, dabs(1.d0-aa/(ONREP*ONREP))
            print*, 'error during the generation of the F_P function is highly &
            &        probable :'
            print*, 'slope of F_P at high freq. - theoretical:', ONREP*OSQNREP
            print*, 'slope of F_P at high freq. - calculated:', aa/OSQNREP
         endif

         ! **** write solution ****
         do j = 1, n
            fP(j) = fP1(j)
            write(100,*) j, x(j)*OSQNREP, fP(j)
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
         dx = dsqrt(dfloat(nrep)) * 0.5d0 * HBokT * dw ! x = u*sqrt(P)
         if (dx > 1.d-3) then
            xmin = 1.d-3
         else
            xmin = dx
         endif
         xmax = dsqrt((n*dx)**2 + dfloat(nrep)**2) ! sqrt(P)*x_kmax
         if (xmax < 50.d0) xmax = 50.d0    ! xmax >= 50 to avoid errors
         nx = int((xmax - xmin)/dx) + 1
         nx = nx + nx/5  !add 20% points to avoid problem at the end
                         !of the interval
         OSQNREP = 1.d0 / dsqrt(dfloat(nrep))
         ONREP1 = 1.d0 / dsqrt(dfloat(nrep) - 1.d0)
         malpha = ONREP ! alpha = 1/P
         niter = 30     ! 30 iterations are enough to converge

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
            call linReg(fP1(9*n/10:n),x(9*n/10:n),aa,bb,r2)

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
         ! "average error"
         if (err > 1.d-6) then
            print*, 'error during the generation of the F_P function is highly &
            &        probable :'
            print*, 'average "error"', err
         endif
         ! value at zero freq.
         if (dabs(1.d0-fP1(1)/ONREP) > 1e-6) then
            print*, 'error during the generation of the F_P function is highly &
            &        probable :'
            print*, 'F_P at zero freq. - theoretical (= 1/P):', ONREP
            print*, 'F_P at zero freq. - calculated:', fP1(1)
         endif
         ! slope at high freq.
         if (dabs(1.d0-aa/(ONREP1*ONREP)) > 5e-3) then
            print*, dabs(1.d0-aa/(ONREP1*ONREP))
            print*, 'error during the generation of the F_P function is highly &
            &        probable :'
            print*, 'slope of F_P at high freq. - theoretical:', ONREP1*OSQNREP
            print*, 'slope of F_P at high freq. - calculated:', aa/OSQNREP
         endif

         ! **** write solution ****
         do j = 1, n
            tmp = dfloat(j) * dx * OSQNREP
            tmp = dfloat(nrep) * (tmp**2 - dfloat(nrep) * dsin(PI*ONREP)**2)
            if (tmp <= 0.d0) then
               fP(j) = 0.0d0
            else
               tmp = dsqrt(tmp)
               i = nint((tmp - xmin) / dx)
               if(i <= 0) then
                  fP(j) = fP1(1)
               else if (i >= n) then
                  fP(j) = fP1(n)
               else
                  fP(j) = fP1(i) + (fP1(i+1)-fP1(i)) * (tmp - x(i))/dx
               endif
            endif
            write(100,*) j, dfloat(j) * dx * OSQNREP, fP(j)
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

      endif

   end subroutine computefP1

   subroutine FFT(x_in, x_out)
      implicit none

      complex*16, dimension(:), intent(in)  :: x_in
      complex*16, dimension(:), intent(out) :: x_out

      integer*8 :: plan ! required for FFTW
      integer :: n, ntest

      n = size(x_in)
      ntest = size(x_out)
      if (n /= ntest) stop('Error with FFT during QTB random force generation')

      ! Forward FFT using FFTW3 library
      call dfftw_plan_dft_1d(plan,n,x_in,x_out,'FFTW_FORWARD','FTTW_ESTIMATE')
      call dfftw_execute_dft(plan,x_in,x_out)
      call dfftw_destroy_plan(plan)

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
         stop('error in linReg subroutine: x and y don''t have the same size')
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

      if (r2 < 0.9d0) print*, 'error in linear regression part of the &
      &f_P function generation is highly probable : r^2 =', r2

   end subroutine linReg

end module thermostats
