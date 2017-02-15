! *************************************************
! simpleMD 1D -  Main program
! *************************************************
! simple MD code for 1 dof systems
! 1 particle in various type of external potentials
! *************************************************
! F. Brieuc - January 2017
! *************************************************
! All units are atomic units
! *************************************************

! contains parameters, constants definition & conversion functions
include 'parameters.f90'
! contains subroutines & functions for I/O
include 'inputOutput.f90'
! contains subroutines & functions for the different thermostats
include 'thermostats.f90'
! contains subroutines & functions to compute forces
include 'forces.f90'
! contains subroutines & functions to compute energies
include 'energies.f90'
! contains subroutines & functions for (pseudo)random numbers generation
include 'ZBQ.f90'

program simpleMD_1d
   use param
   use io
   use force
   use energy
   use thermostats
   implicit none

   character(len = 30) :: inputFile = 'input.dat'
   integer :: istep, irep
   real(dp), dimension(:, :), allocatable :: x  ! position at time t and t-dt
   real(dp), dimension(:), allocatable    :: xx ! "new" position at time t+dt
   real(dp), dimension(:), allocatable    :: v  ! velocitiy at time t
   real(dp), dimension(:), allocatable    :: vv ! "new" velocity at time t+dt
   real(dp) :: F                                ! force that apply on replica i
   real(dp) :: DT2OM, O2DT                      ! dt**2 / (2 * m), 1 / (2 * dt)

   call readInputFile(inputFile)
   call printInfos('b')

   ! usefull values
   DT2OM = dt**2 / m
   O2DT = 1.d0 / (2.d0 * dt)

   allocate(x(nrep, 2))
   allocate(xx(nrep))
   allocate(v(nrep))
   allocate(vv(nrep))

   ! initialization of positions and velocities
   xc = 0.d0
   do irep = 1, nrep
      x(irep, 2) = x0
      xc = xc + x(irep, 2) * ONREP
      v(irep) = v0
   enddo

   ! initialization of the random generator seed
   call ZBQLINI(0)

   ! initialization of thermostats if needed
   call thermostating('initialize')

   ! Taylor step
   call computeInstEner(x(:,2), v)
   xc = 0.d0
   do irep = 1, nrep
      F = totalForce(x(irep,2), x(im(irep),2), x(ip(irep),2), v(irep), irep, istep)
      x(irep,1) = x(irep,2) + dt * v(irep) + (0.5d0 * dt**2 / m) * F
      xc = xc + x(irep, 1) * ONREP
      v(irep) = (x(irep,1) - x(irep,2)) / dt
   enddo

   ! equilibration MD loop
   do istep = 1, neq
      call computeInstEner(x(:,1), v)
      xc = 0.d0
      do irep = 1 ,nrep
         ! forces
         F = totalForce(x(irep,1), x(im(irep),1), x(ip(irep),1), v(irep), irep, istep)
         ! Verlet - position - t+dt
         xx(irep) = 2.d0 * x(irep,1) - x(irep,2) + DT2OM * F
         ! centroid position
         xc = xc + xx(irep) * ONREP
         ! Verlet - velocity - t+dt
         vv(irep) = (3.d0 * xx(irep) - 4.d0 * x(irep,1) + x(irep,2)) * O2DT
      enddo
      ! actualisation t -> t+dt
      do irep = 1, nrep
         x(irep,2) = x(irep,1)
         x(irep,1) = xx(irep)
         v(irep) = vv(irep)
      enddo
   enddo

   ! initialize average energies to zero
   call initAvEner()
   ! averages computation MD loop
   do istep = 1, nstep
      call computeInstEner(x(:,1), v)
      call computeAvEner()
      call printResults(istep, x(:,1))
      xc = 0.d0
      do irep = 1 ,nrep
         ! Force
         F = totalForce(x(irep,1), x(im(irep),1), x(ip(irep),1), v(irep), irep, istep)
         ! Verlet - position - t+dt
         xx(irep) = 2.d0 * x(irep,1) - x(irep,2) + DT2OM * F
         ! centroid position
         xc = xc + xx(irep) * ONREP
         ! Verlet - velocity - t+dt
         vv(irep) = (3.d0 * xx(irep) - 4.d0 * x(irep,1) + x(irep,2)) * O2DT
      enddo
      ! Actualisation t -> t+dt
      do irep = 1, nrep
         x(irep,2) = x(irep,1)
         x(irep,1) = xx(irep)
         v(irep) = vv(irep)
      enddo
   enddo

   deallocate(x)
   deallocate(xx)
   deallocate(v)
   deallocate(vv)

   call printInfos('e')

   ! destroy/close several needed things associated with thermostats if needed
   call thermostating('terminate')
end program simpleMD_1d
