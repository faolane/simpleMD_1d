! *************************************************
! simpleMD 1D -  Module file
! *************************************************
! Contains subroutine and functions to compute forces
! *************************************************
! F. Brieuc - January 2017
! *************************************************

module force
   use param
   implicit none
   private
   public totalForce, potentialForce

contains

   real(dp) function totalForce(x, xm, xp, v, n)
      ! compute the total force that applies on each replica
      implicit none

      real(dp), intent(in) :: x  ! position of replica i
      real(dp), intent(in) :: xm ! position of replica i-1
      real(dp), intent(in) :: xp ! position of replica i+1
      real(dp), intent(in) :: v  ! velocity
      integer, intent(in)  :: n  ! nstep number

      if (nrep > 1) then
         totalForce = potentialForce(x) + thermostatForce(v, n) &
                     &+ replicasForce(x,xm,xp)
      else
         totalForce = potentialForce(x) + thermostatForce(v, n)
      endif

   end function totalForce

   real(dp) function potentialForce(x)
      use param
      implicit none

      real(dp), intent(in) :: x

      select case (pot)
         case('harmonic')
            potentialForce = - m * omega0**2 * x * ONREP
         case('quartic')
         case ('morse')
         case('double-well')
      end select

   end function potentialForce

   real(dp) function thermostatForce(v, n)
      ! compute forces associated with the thermostat
      use param
      use thermostats
      implicit none

      real(dp), intent(in) :: v ! velocity
      integer, intent(in)  :: n ! step number

      select case (therm)
         case('nve')
            thermostatForce = 0.d0
         case('langevin')
            thermostatForce = langevin(v)
         case ('bussi')
         case('qtb')
            thermostatForce = qtb(v, n)
      end select

   end function thermostatForce

   real(dp) function replicasForce(x,xm,xp)
      ! compute forces associated with the PIMD springs btw adjacent replicas
      use param
      implicit none

      real(dp), intent(in) :: x  ! position of replica i
      real(dp), intent(in) :: xm ! position of replica i-1
      real(dp), intent(in) :: xp ! position of replica i+1

      replicasForce = - m * omegaP**2 * (2.d0 * x - xm -xp)

   end function replicasForce

end module force
