! *************************************************
! simpleMD 1D -  Module file
! *************************************************
! Contains subroutines and functions for the
! different thermostats
! *************************************************
! F. Brieuc - January 2017
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
            ! compute random force trough the convolution product
            R = 0.d0
            do i = 0, nQTB-2
               R = R + HQTB(i) * rQTB(irep, i)
               rQTB(irep, i) = rQTB(irep, i + 1)
            enddo
            R = R + HQTB (nQTB - 1) * rQTB(irep, nQTB - 1)
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
      call computefP(fP, dw)
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
      call computefP(fP, dw)
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
      call computefP(fP, dw)
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
      !FFT of H(w) -> H(t)
      call FFT(Hw, Ht)
      do iw = 0, nmax - 1
         HQTB(iw) = realpart(Ht(iw)) * std / dfloat(nmax)
         do idof = 1, nrep
            rQTB(idof, iw) = gaussianRand(1.d0,0.d0)
         enddo
      enddo

      deallocate(Hw)
      deallocate(Ht)

   end subroutine generQTBRandF3

   subroutine computefP(fP, dw)
      use param
      implicit none

      real(dp), intent(in) :: dw
      real(dp), dimension(:), intent(out) :: fP
      integer :: iw, n
      real(dp) :: w, HBokT, tmp

      n = size(fP)
      HBokT = HB / (KB * T)

      if (nrep == 1) then
         do iw = 1, n
            w = dfloat(iw) * dw
            tmp = HBokT * w
            ! fP = theta(w, T) / kT
            fP(iw) = tmp * (0.5d0 + 1.d0 / (exp(tmp) - 1.d0))
         enddo
      endif
   end subroutine computefP

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

end module thermostats
