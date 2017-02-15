! *************************************************
! simpleMD 1D -  Module file
! *************************************************
! Contains parameters and constants definition
! along with other small functions such as for
! unit conversion for example
! *************************************************
! all units are atomic units
! *************************************************
! F. Brieuc - January 2017
! *************************************************

module param
   implicit none

   integer, parameter:: dp = kind(0.d0) ! double precision floating numbers

   ! constants
   real(dp), parameter :: PI = acos(-1.d0)
   real(dp), parameter :: TWOPI = 2.d0 * PI
   real(dp), parameter :: KB = 8.617d-5 / 27.211d0 ! Boltzmann constant (au)
   real(dp), parameter :: HB = 1.d0  ! reduced Planck constant (au)

   ! calculation parameters
   ! thermostats
   character(len = 10) :: therm = 'nve' ! type of thermostat
   real(dp) :: gam = 1.d0               ! friction coeff. in THz (1 THz)
   real(dp) :: T = 300.d0               ! temperature in K
   real(dp) :: fmax = 1.0d7             ! max. freq. (cut-off) for QTB random
                                        ! force generation in THz (1e8 THz!)
   integer :: noiseGenMode = 1          ! noise generation mode
   ! 1 : Random forces are generated using the method of Dammak et al.
   ! 2 : Random forces are generated using a modified version of the method
   ! of Dammak et al. where the random force is generated only for freq. btw
   ! O and fmax and thus the random fore applies for more than one timestep.
   ! 3 : Random forces are generated using the method of Barrat
   ! potential
   integer  :: nQTB = 100   ! nb of points/grid used for the PSD in freq. space
         ! for QTB randf generation with the method of Barrat (noiseGenMode = 3)
   character(len = 15) :: pot = 'harmonic' ! type of external potential
   real(dp) :: f0 = 10.d0                  ! characteristic frequency of
                                           ! harmonic potential in THz (10 THz)
   ! calculation
   integer  :: nstep = 1     ! total number of MD steps (for average values)
   integer  :: neq = 1       ! number of equilibrium steps
   integer  :: nw = 1        ! writing frequency (in outputs files)
   integer  :: nrep = 1      ! number of replicas (Trotter number) for PIMD
   real(dp) :: dt = 1.d0     ! timestep (fs)
   real(dp) :: x0 = 0.d0     ! initial position (bohr)
   real(dp) :: v0 = 0.d0     ! initial velocity (bohr)
   real(dp) :: m = 1.d0      ! mass of the particle (atomic units)
   !other
   integer  :: mQTB = 1      ! rand force generation/update every mQTB steps
   real(dp), dimension(:), allocatable   :: HQTB ! contains the "filter"(Barrat)
   real(dp), dimension(:,:), allocatable :: rQTB ! contains random nb (Barrat)
   real(dp) :: std = 0.d0    ! std deviation of Random force (Lang or QTB)
   real(dp) :: R = 0.d0      ! stores random force (Langevin or QTB)
   real(dp) :: omegamax = 1.d0 ! cut-off ang. freq. for QTB random forces
   real(dp) :: omega0 = 1.d0 ! characteristic ang. freq. of harmonic potential
   real(dp) :: omegaP = 1.d0 ! characteristic ang. freq. of PIMD springs
   real(dp) :: KBTO2 = 1.d0  ! kT/2
   real(dp) :: xc = 0.d0     ! centroid position
   real(dp) :: ONREP = 1.d0  ! 1/nrep
   ! energies
   real(dp) :: Ek = 0.0d0, EkPrim = 0.0d0, EkVir = 0.0d0, EkmVir = 0.0d0
   real(dp) :: Ek_av = 0.0d0, EkPrim_av = 0.0d0, EkVir_av = 0.0d0
   real(dp) :: EkmVir_av = 0.0d0
   ! instantaneous and average kinetic energy computed from different estimators
   ! Ek : "classical", EkPrim : Primitive, EkVir : (centroid) Virial
   ! EkmVir : modified (centroid) Virial. For more infos on these estimators
   ! see Brieuc, Dammak and Hayoun, JCTC, 12, 1351-1359 (2016)
   real(dp) :: Ep = 0.0d0, Ep_av = 0.0d0  ! inst. and av. potential energy

   ! for i/o
   integer :: posFileUnit  = 11
   integer :: velFileUnit  = 12
   integer :: enerFileUnit = 13
   integer :: randfFileUnit = 14

   namelist / thermostat / therm, gam, T, fmax, noiseGenMode, nQTB
   namelist / potential / pot, f0
   namelist / calculation / nstep, neq, nw, nrep, dt, x0, v0, m

contains

   ! converts Hartree to eV
   real(dp) function Ha2eV(ener)
      implicit none

      real(dp), intent(in) :: ener

      Ha2eV = ener * 27.211d0
   end function Ha2eV

   ! converts eV to Hartree
   real(dp) function eV2Ha(ener)
      implicit none

      real(dp), intent(in) :: ener

      eV2Ha = ener / 27.211d0
   end function eV2Ha

   ! converts THZ to atomic frequency unit
   real(dp) function THz2au(freq)
      implicit none

      real(dp), intent(in) :: freq

      THz2au = freq / 4.1341d4
   end function THz2au

   ! converts atomic frequency unit to THZ
   real(dp) function au2THz(freq)
      implicit none

      real(dp), intent(in) :: freq

      au2THz = freq * 4.1341d4
   end function au2THz

   ! converts fs to atomic units
   real(dp) function fs2au(time)
      implicit none

      real(dp), intent(in) :: time

      fs2au = time / 2.4189d-2
   end function fs2au

   ! converts atomic units to fs
   real(dp) function au2fs(time)
      implicit none

      real(dp), intent(in) :: time

      au2fs = time * 2.4189d-2
   end function au2fs

   real(dp) function au2ps(time)
      implicit none

      real(dp), intent(in) :: time

      au2ps = time * 2.4189d-5
   end function au2ps

   integer function ip(i)
      implicit none

      integer, intent(in) :: i

      if (i == nrep) then
         ip = 1
      else
         ip = i + 1
      endif

   end function ip

   integer function im(i)
      implicit none

      integer, intent(in) :: i

      if (i == 1) then
         im = nrep
      else
         im = i - 1
      endif

   end function im

end module param
