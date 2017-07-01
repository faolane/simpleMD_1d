! *************************************************
! simpleMD 1D -  Module file
! *************************************************
! Contains subroutines & functions for the analysis
! of the MD trajectory
! *************************************************
! Copyright (C) 2017 Fabien Brieuc under the terms
! of the GNU General Public License.
! *************************************************

module analysis
   use param
   implicit none
   private
   public analyseMD

contains

   subroutine analyseMD(istep,x)
      use param
      implicit none

      integer, intent(in) :: istep
      real(dp), dimension(nrep), intent(in) :: x  !position

      if (istep == 1) then
         if (boolProba) call computeProbaDens(x,0)
      endif

      if(boolProba) call computeProbaDens(x,1)
      if(boolRadGyr) call computeRadiusofGyration(x)

      if (istep == nstep) then
         if (boolProba) call computeProbaDens(x,2)
      endif
   end subroutine analyseMD

   subroutine computeProbaDens(x,mode)
      ! compute the proba. density of finding the particle at position x
      use param
      implicit none

      real(dp), dimension(nrep), intent(in) :: x  !position
      integer,  intent(in) :: mode   ! mode = 0 (init), 1(running), 2(end)
      real(dp) :: norm ! normalization
      integer :: i,k

      if (mode == 0) then
         allocate(proba(0:Ndens))
         allocate(probaCen(0:Ndens))
         proba = 0.d0
         probaCen = 0.d0
         dxdens = (x1dens - x0dens) / Ndens
      else if (mode == 1) then
         do i = 1, nrep
            k = nint((x(i) - x0dens) / dxdens)
            if (k >= 0 .and. k <= Ndens) proba(k) = proba(k) + 1
         enddo
         k = nint((xc - x0dens) / dxdens)
         if (k >= 0 .and. k <= Ndens) probaCen(k) = probaCen(k) + 1
      else if (mode == 2) then
         ! normalization
         norm = sum(proba) * dxdens
         proba = proba / norm
         norm = sum(probaCen) * dxdens
         probaCen = probaCen / norm
         do i = 0, Ndens
            write(probaFileUnit,*) x0dens + i * dxdens, proba(i), probaCen(i)
         enddo
         deallocate(proba)
         deallocate(probaCen)
      endif

   end subroutine computeProbaDens

   subroutine computeRadiusofGyration(x)
      ! compute the radius of gyration
      use param
      implicit none

      real(dp), dimension(nrep), intent(in) :: x  !position
      real(dp) :: tmp
      integer :: i

      tmp = 0.d0
      do i = 1, nrep
         tmp = tmp + (x(i) - xc)**2
      enddo
      tmp = tmp*ONREP
      radGyr = radGyr + dsqrt(tmp)

   end subroutine computeRadiusofGyration

end module analysis
