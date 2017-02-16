! *************************************************
! simpleMD 1D -  Module file
! *************************************************
! Contains subroutine and functions to compute energies
! *************************************************
! Copyright (C) 2017 Fabien Brieuc under the terms
! of the GNU General Public License.
! *************************************************

module energy
   use param
   implicit none
   private
   public initAvEner, computeInstEner, computeAvEner

contains

   subroutine initAvEner()
      use param
      implicit none

      Ek_av = 0.d0
      EkPrim_av = 0.d0
      EkVir_av = 0.d0
      EkmVir_av = 0.d0
      Ep_av = 0.d0

   end subroutine initAvEner

   subroutine computeInstEner(x, v)
      implicit none

      real(dp), dimension(nrep), intent(in) :: x  !position
      real(dp), dimension(nrep), intent(in) :: v  !velocity
      integer :: i
      real(dp) :: Ek_tmp, EkPrim_tmp, EkVir_tmp, EkmVir_tmp, Ep_tmp
      real(dp) :: vc

      Ek_tmp = 0.d0; EkPrim_tmp = 0.d0; EkVir_tmp = 0.d0; EkmVir_tmp = 0.d0
      Ep_tmp = 0.d0

      vc = 0.d0
      do i = 1, nrep
         vc = vc + v(i) * ONREP
         Ek_tmp = Ek_tmp + Ekin(v(i))
         EkPrim_tmp = EkPrim_tmp + EkinPrim(x(i), x(ip(i)))
         EkVir_tmp = EkVir_tmp + EkinVir(x(i))
         Ep_tmp = Ep_tmp + Epot(x(i))
      enddo

      Ek = Ek_tmp
      if (trim(therm) == 'qtb') then
         EkPrim = Ek_tmp + EkPrim_tmp
         EkVir = Ek_tmp * ONREP + EkVir_tmp
         EkmVir = nrep * Ekin(vc) + EkVir_tmp
      else
         EkPrim = nrep * KBTO2 + EkPrim_tmp
         EkVir = KBTO2 + EkVir_tmp
         EkmVir = nrep * Ekin(vc) + EkVir_tmp
      endif
      Ep = Ep_tmp

   end subroutine computeInstEner

   subroutine computeAvEner()
      use param
      implicit none

      Ek_av = Ek_av + Ek
      EkPrim_av = EkPrim_av + EkPrim
      EkVir_av = EkVir_av + EkVir
      EkmVir_av = EkmVir_av + EkmVir
      Ep_av = Ep_av + Ep

   end subroutine computeAvEner

   real(dp) function Ekin(v)
      use param
      implicit none

      real(dp), intent(in) :: v ! velocity

      Ekin = 0.5d0 * m * v**2

   end function Ekin

   real(dp) function EkinPrim(x, xp)
      use param
      implicit none

      real(dp), intent(in) :: x ! position
      real(dp), intent(in) :: xp ! position

      EkinPrim = - 0.5d0 * m * omegaP**2 * (x - xp)**2

   end function EkinPrim

   real(dp) function EkinVir(x)
      use param
      use force
      implicit none

      real(dp), intent(in) :: x ! position

      EkinVir = - (x - xc) * potentialForce(x) * 0.5D0

   end function EkinVir

   real(dp) function Epot(x)
      use param
      implicit none

      real(dp), intent(in) :: x ! position

      select case (pot)
         case('harmonic')
            Epot = 0.5d0 * m * omega0**2 * x**2 * ONREP
         case('quartic')
            Epot = 0.0d0
         case ('morse')
            Epot = 0.0d0
         case('double-well')
            Epot = 0.0d0
         case default
            Epot = 0.0d0
      end select
   end function Epot

end module energy
