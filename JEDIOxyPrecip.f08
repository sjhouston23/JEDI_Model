subroutine JEDIOxyPrecip(version,Oxy)
!*******************************************************************************
!* Created by Stephen J. Houston 9.4.18
!*******************************************************************************
!* This program reads in a JEDI spectrum (*.d2s) and normalizes it based on the
!* energy and the energy bin width. It then multiplies the new input ion flux
!* by normalized (1 input ion/cm^2/s), previously calculated results that
!* correspond to the JEDI energy bins ~(11,15,20,30,47,60,78,121,218,456 keV).
!* Some of the higher energy JEDI bins overlap - I treat them as if they don't.
!* To change the number of interpolated data points, one would need to change
!* the variables number_of_energies, the Eion data points, and then go into
!* the JEDIInterpolator subroutine and make changes.
!*******************************************************************************

use, intrinsic :: ISO_FORTRAN_ENV !Used for kind=int64
use formattingOxy !Used for formatting.f08
implicit real*8(a-h,o-z) !i,j,k,l,m,n are assumed to be integers

!*****************INTERFACE BLOCK FOR ALLOCATABLE ARRAYS*****************
!You must add an interface block for any allocatable arrays that you want to
!pass from one routine/subroutine to another.

interface
  subroutine FluxInterpolator(PJ,spec,Ienergy,Iflux)
    character(len=100),intent(in) :: PJ
    integer,intent(in) :: spec
    real*8,allocatable,dimension(:),intent(out) :: Ienergy,Iflux

    integer oxy,sul
    parameter(nSteps=3,oxy=1,sul=2)
    !nSteps is the number of energy data points you want in between the JEDI points

    !* JEDI variables
    real*8 JebinsMIN,JebinsMAX
    parameter(JebinsMIN=135.62,JebinsMAX=10000.00) !From Mauk et al. (2013) - SSR
    real*8,allocatable,dimension(:) :: Jenergy,Jintensity,Jebins,Jflux,JebinWidth
    !*   Jenergy - JEDI energy bin [keV]
    !*   Jintensity - JEDI ion flux [c/s/ster/cm^2/keV]
    !*   Jebins - Size of JEDI energy bins [keV]
    !*   Jflux - Jintensity values converted to [counts/cm^2/s]

    ! character(len=10) date
    ! character(len=12) time
    character(len=100) filename,version
    character*100 tmp
    !* Interpolated variables
    real*8,allocatable,dimension(:) :: Iintensity,Iebins
  end subroutine
end interface

!*******************************************************************************
character(len=100) :: version
real*8,intent(out) :: Oxy !Total flux of oxygen ions

integer energy,atmosLen,run

parameter(atmosLen=1544,nProc=36,nChS=10) !Atmosphere, processes, charge states
parameter(nE2strBins=260) !Number of 2-stream energy bins
parameter(nOutputFiles=15) !Number of data files from ion precip code
parameter(nEnergyFiles=53) !Number of energies from ion precip code
! parameter(number_of_energies=39) !Number of interpolated JEDI energy bins
integer,dimension(nEnergyFiles) :: Eprecip !Precipitation code energies
real*8 diff(nEnergyFiles) !Used to find which energy files should be used
! real*8 Eion(number_of_energies) !Ion energies
real*8,dimension(atmosLen) :: altitude !Altitude array
real*8,allocatable,dimension(:,:) :: totalHp,totalH2p,H2Ex
real*8,allocatable,dimension(:,:,:) :: Hp,H2p
real*8,allocatable,dimension(:,:,:,:) :: oxygen
real*8,allocatable,dimension(:,:,:) :: prode2stF,prode2stB
!* JEDI variables
real*8,allocatable,dimension(:) :: Jflux,Eion
!*   Jflux - JEDI measured flux intensities converted to [counts/cm^2/s]

!Same as previous variables except has a leading "J"
real*8 IntegratedPhot(2,nChS) !1 is DE, 2 is CX, integrated photon production
real*8,dimension(atmosLen) :: JtotalHp,JtotalH2p,JH2Ex
real*8,dimension(atmosLen,nProc) :: JHp,JH2p
real*8,dimension(atmosLen,nProc,nChS) :: Joxygen
real*8,dimension(atmosLen,nE2strBins) :: Jprode2stF,Jprode2stB

character(len=100) filename,filenames(nOutputFiles)
character(len=1000) HpHeader,Hp2Header

!****************************** Data Declaration *******************************
!* Initial ion enegy input:
!data Eion/10.625,15.017,20.225,29.783,46.653,59.770,77.522,120.647,218.125,&
!          456.250/ !Original Juno energy bins from JEDI.
!* Energies from the precipitation code (keV/u)
data Eprecip/10,11,12,13,14,15,16,17,19,20,22,25,27,30,33,37,42,47,50,53,56,60,&
             64,68,73,75,78,87,97,100,108,121,125,140,150,162,175,188,200,218,&
             250,262,300,315,350,379,400,450,456,500,541,600,625/
! data Eion/10.625,11.619,12.656,13.786,15.017,16.177,17.427,18.774,20.225,&
! 22.280,24.543,27.036,29.783,33.319,37.276,41.702,46.653,49.634,52.806,56.180,&
! 59.770,63.785,68.070,72.642,77.522,86.586,96.710,108.018,120.647,139.899,&
! 162.223,188.108,218.125,262.319,315.467,379.384,456.250,541.463,625/ !JEDI
data filenames/'H+_Prod','H2+_Prod','H2_Excite_Prod','Oxy_Neg','Oxy0_','Oxy1_',&
'Oxy2_','Oxy3_','Oxy4_','Oxy5_','Oxy6_','Oxy7_','Oxy8_','2Str_Elect_Fwd',&
'2Str_Elect_Bwd'/ !Filenames that are read in from JupOxyPrecip code
!* Width of JEDI energy bins: (May eventually need to be adjusted)
!data Jebins/66.0,71.0,105.0,216.0,346.0,251.0,300.0,880.0,2280.0,5340.0/
!********************************* Initialize **********************************
pi=4.0d0*atan(1.0d0)
!*************************** Call JEDI Interpolator ****************************
! version='PJ7_Paper'
! write(version,'("v1")') !Filename of a JEDI spectrum (.d2s file)
call FluxInterpolator(version,1,Eion,Jflux)
Oxy=sum(Jflux) !Used to calculat O:S ratio later
! goto 10001
!*** Allocation of variables
allocate(totalHp(size(Eion),atmosLen),totalH2p(size(Eion),atmosLen))
allocate(H2Ex(size(Eion),atmosLen))
allocate(Hp(size(Eion),atmosLen,nProc),H2p(size(Eion),atmosLen,nProc))
allocate(oxygen(size(Eion),atmosLen,nProc,nChS))
allocate(prode2stF(size(Eion),atmosLen,nE2strBins))
allocate(prode2stB(size(Eion),atmosLen,nE2strBins))
!********************************* Initialize **********************************
energy=0;altitude=0.0;Hp=0.0;totalHp=0.0;H2p=0.0;totalH2p=0.0;H2Ex=0.0
oxygen=0.0;prode2stF=0.0;prode2stB=0.0
!********** Open output data files for each set of initial energies ************
write(*,*) 'Opening oxygen preciptation files for all energies...'
do run=1,size(Eion) !Loop through each initial ion energy
  !Find the closest precipitation code energy
  do i=1,nEnergyFiles
    diff(i)=abs((Eion(run)/16)-Eprecip(i)) !Find the difference in energy
  end do
  energy=Eprecip(minloc(diff,1)) !Use the minimum difference location (index)
  do i=1,nOutputFiles !Open all of the files
    write(filename,'("./PrecipFiles/Oxy/",I0,"keV/",A,"_Comb.dat")') &
          energy,trim(filenames(i))
    filename=trim(filename)
    open(unit=100+i,file=filename,status='old')
    rewind(100+i)
  end do
  do i=1,2
    read(101,*) !Read the headers of the first four files
    read(102,*)
    read(103,*)
    read(103,*)
  end do
  read(101,'(1x,A)') HpHeader !Save these headers to output later
  read(102,'(1x,A)') Hp2Header
  do i=1,atmosLen !Hydrogen Ionization/Excitation vs. altitude
    read(101,F01) altitude(i),(Hp(run,i,j),j=1,31),totalHp(run,i)
    read(102,F01) dum,(H2p(run,i,j),j=1,11),totalH2p(run,i)
    read(103,F02) dum,H2Ex(run,i)
  end do
  do i=1,nChS !Loop through every charge state
    read(103+i,*) !Oxygen header
    do j=1,atmosLen
      read(103+i,F01) dum,(oxygen(run,j,k,i),k=1,nProc)
    end do
  end do
  do j=1,nE2strBins !Electron production rate for 2-stream
    read(114,F2Str) (prode2stF(run,i,j),i=atmosLen,1,-1)
    read(115,F2str) (prode2stB(run,i,j),i=atmosLen,1,-1)
  end do
end do !End of do-loop for each energy
do i=1,nOutputFiles
  close(100+i) !close files
end do
!********************************* Initialize **********************************
JHp=0.0;JtotalHp=0.0;JH2p=0.0;JtotalH2p=0.0;JH2Ex=0.0;Joxygen=0.0
Jprode2stF=0.0;Jprode2stB=0.0;IntegratedPhot=0.0
!************************* Calculate JEDI Productions **************************
write(*,*) 'Calculating JEDI production rates...'
do run=1,size(Eion) !Loop through every energy bin
  do i=1,atmosLen !Loop through entire atmosphere
    JtotalHp(i)=JtotalHp(i)+totalHp(run,i)*Jflux(run) !Total H+
    JtotalH2p(i)=JtotalH2p(i)+totalH2p(run,i)*Jflux(run) !Total H2+
    JH2Ex(i)=JH2Ex(i)+H2Ex(run,i)*Jflux(run) !Total H2*
    do j=1,nProc !Loop through every processes
      JHp(i,j)=JHp(i,j)+Hp(run,i,j)*Jflux(run) !H+ by process
      JH2p(i,j)=JH2p(i,j)+H2p(run,i,j)*Jflux(run) !H2+ by process
      do k=1,nChS !Loop through every charge state
        Joxygen(i,j,k)=Joxygen(i,j,k)+Oxygen(run,i,j,k)*Jflux(run) !Oxygen
      end do !End charge state loop - k
    end do !End processes loop - j
    do j=1,nE2strBins !Loop through 2-stream energy bins
      Jprode2stF(i,j)=Jprode2stF(i,j)+prode2stF(run,i,j)*Jflux(run) !e- forward
      Jprode2stB(i,j)=Jprode2stB(i,j)+prode2stB(run,i,j)*Jflux(run) !e- backward
    end do !End 2-stream energy bins loop
  end do !End atmsophere loop
end do !End energy bin loop
!************************* Write Out JEDI Productions **************************
write(*,*) 'Writing output files...'
do i=1,nOutputFiles
  write(filename,'("./Output/Oxy/",A,"/",A,".dat")') &
        trim(version),trim(filenames(i))
  filename=trim(filename)
  open(unit=200+i,file=filename,status='unknown')
end do
write(201,H01) !H+ header
write(201,*) trim(HpHeader)
write(202,H02) !H2+ header
write(202,*) trim(Hp2Header)
write(203,H03) !H2* header
do i=1,atmosLen !Ionization/Excitation vs. altitude
  write(201,F01) altitude(i),(JHp(i,j),j=1,31),JtotalHp(i)
  write(202,F01) altitude(i),(JH2p(i,j),j=1,11),JtotalH2p(i)
  write(203,F02) altitude(i),JH2Ex(i)
end do
do i=1,nChS !Oxygen production
  write(203+i,*) "Alt [km] ", (HProc(k),k=1,nProc)
  do j=1,atmosLen
    write(203+i,F01) altitude(j),(Joxygen(j,k,i),k=1,nProc)
  end do
end do
!****************************** Photon Production ******************************
write(filename,'("./Output/Oxy/",A,"/XRay_DE.dat")') trim(version)
open(unit=301,file=trim(filename),status='unknown') !Open photon DE
write(filename,'("./Output/Oxy/",A,"/XRay_CX.dat")') trim(version)
open(unit=302,file=trim(filename),status='unknown') !Open photon CX
write(filename,'("./Output/Oxy/",A,"/Total_Photon_Prod.dat")') trim(version)
open(unit=303,file=trim(filename),status='unknown') !Open total photon prod
altDelta=2.0e5 !2 km = 200,000 cm
write(301,N01) !DE X-Ray note
write(302,N02) !CX X-Ray note
write(301,*) !DE X-Ray note
write(302,*) !CX X-Ray note
do i=1,3
  ! write(300+i,H11) !Initial input
  ! write(300+i,*) !Extra space
  write(300+i,H08) !Altitude integrated photon production header
  write(300+i,H09) !Charge state header
end do
do i=1,atmosLen !DE - TEX+SPEX,SI+SPEX,DI+SPEX, CX - SC+SS,TI,SC
  do j=1,nChS
    IntegratedPhot(1,j)=IntegratedPhot(1,j)+&
      real(Joxygen(i,27,j)+Joxygen(i,29,j)+Joxygen(i,32,j))*altDelta !DE
    IntegratedPhot(2,j)=IntegratedPhot(2,j)+&
      real(Joxygen(i,19,j)+Joxygen(i,25,j)+Joxygen(i,30,j))*altDelta !CX
  end do
end do
write(303,F05) altDelta/1e5,((IntegratedPhot(1,j)+IntegratedPhot(2,j)),j=1,nChS)
do i=1,3
  if(i.le.2)write(300+i,F05) altDelta/1e5,(IntegratedPhot(i,j),j=1,nChS)
  write(300+i,*) !Extra space
  write(300+i,H10) !Photon production vs. altitude header
  write(300+i,H06) !Charge state header
end do
do i=1,atmosLen !DE - TEX+SPEX,SI+SPEX,DI+SPEX, CX - SC+SS,TI,SC
  write(301,F05) altitude(i),& !Photon production from direct excitation
    ((Joxygen(i,27,j)+Joxygen(i,29,j)+Joxygen(i,32,j)),j=1,nChs)
  write(302,F05) altitude(i),& !Photon production from charge exchange
    ((Joxygen(i,19,j)+Joxygen(i,25,j)+Joxygen(i,30,j)),j=1,nChs)
  write(303,F05) altitude(i),& !Total photon production
    ((Joxygen(i,27,j)+Joxygen(i,29,j)+Joxygen(i,32,j)+&
      Joxygen(i,19,j)+Joxygen(i,25,j)+Joxygen(i,30,j)),j=1,nChs)
end do
close(301)
close(302)
close(303)
!***************************** Secondary Electrons *****************************
do j=1,nE2strBins !2-Stream electrons, forward and backward
  write(214,F2Str) (Jprode2stF(i,j),i=atmosLen,1,-1)
  write(215,F2Str) (Jprode2stB(i,j),i=atmosLen,1,-1)
end do
do i=1,nOutputFiles !Close the combine output files
  close(200+i)
end do
! 10001 continue
end subroutine
