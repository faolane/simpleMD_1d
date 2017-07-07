! *************************************************
! simpleMD 1D -  Module file
! *************************************************
! Contains parameters and constants definition
! along with other small functions such as for
! unit conversion for example
! *************************************************
! all units are atomic units
! *************************************************
! Copyright (C) 2017 Fabien Brieuc under the terms
! of the GNU General Public License.
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
   real(dp) :: fcut = 1.0d7             ! max. freq. (cut-off) for QTB random
                                        ! force generation in THz (1e7 THz!)
   integer :: noiseGenMode = 1          ! noise generation mode
   ! 1 : Random forces are generated using the method of Dammak et al.
   ! 2 : Random forces are generated using a modified version of the method
   ! of Dammak et al. where the random force is generated only for freq. btw
   ! O and fcut and thus the random fore applies for more than one timestep.
   ! 3 : Random forces are generated using the method of Barrat
   ! potential
   integer  :: nQTB = 250   ! nb of points/grid used for the PSD in freq. space
         ! for QTB randf generation with the method of Barrat (noiseGenMode = 3)
   integer  :: piqtbMode = 0               ! QTB-PIMD mode (0)
   character(len = 15) :: pot = 'harmonic' ! type of external potential
   real(dp) :: f0 = 10.d0                  ! characteristic frequency of
                                           ! harmonic potential in THz (10 THz)
   real(dp) :: A = 1.d0                    ! for quartic potential V(x)=A*x^4
   real(dp) :: V0 = 1.d0                   ! barrier height of the double-well
                                           ! potential in eV
   real(dp) :: x0 = 1.d0                   ! position of the two wells (in bohr)
   real(dp) :: D = 1.d0                    ! Dissociation energy of the Morse
                                           ! potential (in eV)
   real(dp) :: alpha = 1.d0                ! width of the Morse pot. (1/bohr)
   ! calculation
   integer  :: nstep = 1      ! total number of MD steps (for average values)
   integer  :: neq = 1        ! number of equilibrium steps
   integer  :: nw = 1         ! writing frequency (in outputs files)
   integer  :: nrep = 1       ! number of replicas (Trotter number) for PIMD
   real(dp) :: dt = 1.d0      ! timestep (fs)
   real(dp) :: xini = 0.d0    ! initial position (bohr)
   real(dp) :: vini = 0.d0    ! initial velocity (bohr)
   real(dp) :: m = 1.d0       ! mass of the particle (atomic units)
   logical  :: auto = .false. ! automatic computation of calculation parameters
   real(dp) :: gamOfmin = 0.2d0 ! gamma / freq. min.
   real(dp) :: dtOTmin = 0.02d0 ! dt / periode min.
   real(dp) :: gamNdt = 4000.d0 ! gamma * Nmax * dt
   real(dp) :: fcutOfmax = 2.d0 ! fcut / fmax
   real(dp) :: QTBfreqStep = 0.5d0 ! freq step used for generation of H(w) and
                                   ! H(t) in QTB - noiseGenMode = 3 (THz)

   ! other
   integer  :: cpt = 10      ! to follow the trajectory advancement
   integer  :: mQTB = 1      ! rand force generation/update every mQTB steps
   real(dp), dimension(:), allocatable   :: HQTB ! contains the "filter"(Barrat)
   real(dp), dimension(:,:), allocatable :: rQTB ! contains random nb (Barrat)
   real(dp) :: std = 0.d0    ! std deviation of Random force (Lang or QTB)
   real(dp) :: omegacut = 1.d0 ! cut-off ang. freq. for QTB random forces
   real(dp) :: omega0 = 1.d0 ! characteristic ang. freq. of harmonic potential
   real(dp) :: omegaP = 1.d0 ! characteristic ang. freq. of PIMD springs
   real(dp) :: KBTO2 = 1.d0  ! kT/2
   real(dp) :: MW02          ! m * omega0**2
   real(dp) :: TDALP         ! 2.d0 * D * alpha
   real(dp) :: FV0oX02       ! 4.d0 * V0 / x0**2
   real(dp) :: xc = 0.d0     ! centroid position
   real(dp) :: ONREP = 1.d0  ! 1/nrep
   real(dp),dimension(:),allocatable :: R ! stores random forces (Langevin, QTB)
   real(dp),dimension(:),allocatable :: Rk ! random forces on NM (QTB-PIMD)
   real(dp),dimension(:,:),allocatable :: C ! Transformation matrix from NM
   ! coordinates to direct coordinates for random forces
   ! energies
   real(dp) :: Ek = 0.0d0, EkPrim = 0.0d0, EkVir = 0.0d0, EkmVir = 0.0d0
   real(dp) :: Ec = 0.0d0
   real(dp) :: Ek_av = 0.0d0, EkPrim_av = 0.0d0, EkVir_av = 0.0d0
   real(dp) :: EkmVir_av = 0.0d0, Ec_av = 0.0d0
   ! instantaneous and average kinetic energy computed from different estimators
   ! Ek : "classical", EkPrim : Primitive, EkVir : (centroid) Virial
   ! EkmVir : modified (centroid) Virial. For more infos on these estimators
   ! see Brieuc, Dammak and Hayoun, JCTC, 12, 1351-1359 (2016)
   real(dp) :: Ep = 0.0d0, Ep_av = 0.0d0  ! inst. and av. potential energy

   ! analysis
   ! proba
   logical  :: boolProba = .False.
   real(dp) :: x0dens = 0.0d0 ! start of the interval for proba. dens.
   real(dp) :: x1dens = 1.0d0 ! end of the interval for proba. dens.
   real(dp) :: dxdens         ! dx for computation of proba (width of the bins)
   integer  :: Ndens = 100    ! number of bins
   real(dp), dimension(:), allocatable :: proba ! replicas proba. density
   real(dp), dimension(:), allocatable :: probaCen ! centroid proba. density
   ! radius of gyration
   logical  :: boolRadGyr = .False.
   real(dp) :: radGyr = 0.d0 ! radius of gyration

   ! for i/o
   integer :: posFileUnit  = 11
   integer :: velFileUnit  = 12
   integer :: enerFileUnit = 13
   integer :: randfFileUnit = 14
   integer :: probaFileUnit = 15
   integer :: radGyrFileUnit = 16

   namelist / thermostat / therm, gam, T, fcut, noiseGenMode, nQTB, piqtbMode
   namelist / potential / pot, f0, A, V0, x0, D, alpha
   namelist / calculation / nstep, neq, nw, nrep, dt, xini, vini, m, auto,&
                            gamOfmin, dtOTmin, gamNdt, fcutOfmax, QTBfreqStep
   namelist / analyse / boolProba, x0dens, x1dens, Ndens, boolRadGyr

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
