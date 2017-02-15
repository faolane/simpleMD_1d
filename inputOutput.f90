! *************************************************
! simpleMD 1D -  Module file
! *************************************************
! Contains subroutines for I/O
! *************************************************
! Copyright (C) 2017 Fabien Brieuc under the terms
! of the GNU General Public License.
! *************************************************

module io
   use param
   implicit none
   private
   public readInputFile, printInfos, printResults

contains

   subroutine readInputFile(fileName)
      ! read input parameters, initialize and converts it
      use param
      implicit none

      character(len = 30), intent(in) :: fileName

      open(10, file = trim(fileName))
      read(10, calculation)
      dt = fs2au(dt)
      read(10, potential)
      omega0 = TWOPI * THz2au(f0)
      read(10, thermostat)
      gam = THz2au(gam)
      fmax = THz2au(fmax)
      omegamax = TWOPI * fmax
      close(10)

      ! characteristic angular freq. of the PIMD springs btw replicas
      if (nrep > 1) then
         omegaP = dsqrt(dfloat(nrep)) * KB * T / HB
      else
         omegaP = 0.d0
      endif

      ! some usefull quantities
      KBTO2 = 0.5d0 * KB * T
      ONREP = 1.d0 / dfloat(nrep)

   end subroutine readInputFile

   subroutine printInfos(mode)
      ! print infos on the screen
      use param
      implicit none

      character, intent(in) :: mode
      integer, dimension(3) :: date, time ! store date and time

      select case (mode)
         case('b')
            call idate(date); call itime(time)
            write(6, *) '###### simpleMD 1D ######'
            write(6, '(a25,i2.2,a1,i2.2,a1,i4,a4,i2.2,a,i2.2,a,i2.2)'), &
                     &' calculation launched on ', &
                     & date(1), '\',  date(2), '\', date(3), &
                     &' at ', time(1), ':', time(2), ':', time(3)
            if (nrep > 1) then
               write(6, *) ''
               write(6, *) '# PIMD Simulation #'
               write(6, '(a28,i4)') ' number of replicas : nrep = ', nrep
               write(6, '(a16,e11.2,a3)') ' timestep : dt = ', au2fs(dt), ' fs'
               write(6, '(a42,i8)') ' number of equilibration MD steps : &
                         &neq = ', neq
               write(6, '(a50,i10)') ' number of average computation MD &
                         &steps : nstep = ', nstep
               write(6, '(a19,f7.2,a2)') ' temperature : T = ', T, ' K'
            else
               write(6, *) ''
               write(6, *) '# MD Simulation #'
               write(6, '(a28,i4)') ' number of replicas : nrep = ', nrep
               write(6, '(a16,e11.2,a3)') ' timestep : dt = ', au2fs(dt), ' fs'
               write(6, '(a42,i8)') ' number of equilibration MD steps : &
                         &neq = ', neq
               write(6, '(a50,i10)') ' number of average computation MD &
                         &steps : nstep = ', nstep
               write(6, '(a19,f7.2,a2)') ' temperature : T = ', T, ' K'
            endif
            write(6, *) ''
            select case (therm)
               case('nve')
                  write(6, *) '# NVE (No thermostat) #'
               case('qtb')
                  write(6, *) '# QTB thermostat #'
                  write(6, '(a29,e11.2,a4)') ' friction coefficient : gam = ', &
                                              &au2THz(gam), ' THz'
                  write(6, '(a27,e11.2,a4)') ' cut-off frequency : fmax = ', &
                                              &au2THz(fmax), ' THz'
               case('bussi')
                  write(6, *) '# Bussi thermostat #'
               case('langevin')
                  write(6, *) '# Langevin thermostat #'
                  write(6, '(a29,e11.2,a4)') ' friction coefficient : gam = ', &
                                              &au2THz(gam), ' THz'
            end select
            write(6, *) ''
            select case (pot)
               case('harmonic')
                  write(6, *) '# Harmonic potential #'
                  write(6, '(a28,f7.2,a4)') ' characteristic frequency : ', &
                                             &au2THz(omega0 / TWOPI), ' THz'

               case('morse')
                  write(6, *) '# Morse potential #'
               case('double-well')
                  write(6, *) '# Double-well potential #'
               case('quartic')
                  write(6, *) '# Quartic potential #'
               end select
         case('e')
            print*, ''
            call idate(date); call itime(time)
            write(6,'(a35,i2.2,a1,i2.2,a1,i4,a4,i2.2,a,i2.2,a,i2.2)'), &
                     &' calculation correctly finished on ', &
                     & date(1), '\',  date(2), '\', date(3), &
                     &' at ', time(1), ':', time(2), ':', time(3)
      end select

   end subroutine printInfos

   subroutine printResults(i, x)
      ! print results in the different associated files
      use param
      implicit none

      integer, intent(in) :: i
      real(dp), dimension(nrep), intent(in) :: x ! positions

      ! open files
      if (i == 1) then
         ! positions file
         open(posFileUnit, file='positions.res')
         write(posFileUnit, *) '# step nb, time (ps), position rep 1 (bohr), &
                     &..., position rep nrep (bohr), centroid position (bohr)'
         ! velocities file
         ! open(velFileUnit, file='velocities.res')
         ! write(velFileUnit, *) '# step nb, time(ps), velocity rep 1 (au), &
         !                &..., velocity rep nrep (au)'
         ! energies file
         open(enerFileUnit, file='energies.res')
         if (nrep > 1) then
            write(enerFileUnit, *) '# step nb, time (ps), ''classical'' kinetic &
               &energy (eV), average ''classical'' kin. ener. (eV), primitive &
               &kin. ener. (eV), av. prim. kin. ener. (eV), virial kin. ener. &
               &(eV), av. vir. kin. ener. (eV), modified vir. kin. ener. (eV), &
               &av. mod. vir. kin. ener. (eV), potential ener. (eV), av. pot. &
               &ener. (eV)'
         else
            write(enerFileUnit, *) '# step nb, time (ps), kinetic &
               &energy (eV), average kin. ener. (eV), potential &
               &ener. (eV), av. pot. ener. (eV)'
         endif
      endif

      if (mod(i, nw) == 0) then
         ! print positions
         write(posFileUnit, *) i, au2ps(i*dt), x, xc
         ! print velocities
         ! write(velFileUnit, *) i, i*dt, v
         ! print energies
         if (nrep > 1) then
            write(enerFileUnit, *) i, au2ps(i*dt), Ha2eV(Ek), Ha2eV(Ek_av) / i, &
                                  &Ha2eV(EkPrim), Ha2eV(EkPrim_av) / i, &
                                  &Ha2eV(EkVir), Ha2eV(EkVir_av) / i, &
                                  &Ha2eV(EkmVir), Ha2eV(EkmVir_av) / i, &
                                  &Ha2eV(Ep), Ha2eV(Ep_av) / i
         else
            write(enerFileUnit, *) i, au2ps(i*dt), Ha2eV(Ek), Ha2eV(Ek_av) / i, &
                                  &Ha2eV(Ep), Ha2eV(Ep_av) / i
         endif
      endif

      ! close files
      if (i == nstep) then
         close(posFileUnit)
         ! close(velFileUnit)
         close(enerFileUnit)
      endif

   end subroutine printResults
end module io
