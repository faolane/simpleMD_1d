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
         proba = 0.d0
         dxdens = (x1dens - x0dens) / Ndens
         print*, 'dx', dxdens
      else if (mode == 1) then
         do i = 1, nrep
            k = nint((x(i) - x0dens) / dxdens)
            if (i >=0 .and. i<=Ndens) proba(k) = proba(k) + 1
         enddo
      else if (mode == 2) then
         ! normalization
         norm = sum(proba) * dxdens
         print*, 'norm', norm
         proba = proba / norm
         open(200, file='proba-density.res')
         write(200,*) '# probability density of position (1/bohr)'
         write(200,*) '# position x (bohr), proba density (1/bohr)'
         do i = 0, Ndens
            write(200,*) x0dens + i * dxdens, proba(i)
         enddo
         close(200)
         deallocate(proba)
      endif

   end subroutine computeProbaDens

end module analysis
