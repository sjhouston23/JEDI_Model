! program FluxInterpolator
subroutine FluxInterpolator(PJ,spec,Ienergy,Iflux)
!*******************************************************************************
!* Created by Stephen J. Houston 9.5.18
!*******************************************************************************
!* This program interpolates JEDI ion spectra into a finer energy grid. It
!* outputs the intensity, energy bin width, and then a renormalized ion flux
!* which is too be input into the oxygen ion precipitation model. Alternatively,
!* the flux can instead be multiplied by a normalized (1 ion/cm^2/s) output from
!* JupOxyPrecip.
!*******************************************************************************

implicit real*8(a-h,o-z)

!*******************************************************************************
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

!****************************** Data Declaration *******************************
!* Width of JEDI energy bins: (May eventually need to be adjusted)
! data Jebins/66.0,71.0,105.0,216.0,346.0,251.0,300.0,880.0,2280.0,5340.0/
!********************************* Initialize **********************************
pi=4.0d0*atan(1.0d0)
! spec=sul
! PJ='PJ7_Paper'
! Jenergy=0.0;Jintensity=0.0;Jbins=0.0;Jflux=0.0
!*************************** Open JEDI Ion Spectrum ****************************
if(spec.eq.oxy)then
  write(version,'("Oxy/",A)') trim(PJ)
elseif(spec.eq.sul)then
  write(version,'("Sul/",A)') trim(PJ)
end if
! write(version,'("Oxy/PJ7_Paper")') !Filename of a JEDI spectrum (.d2s file)
write(filename,'("./JunoData/Spectra/",A,".d2s")') trim(version)
write(*,*) filename
open(unit=100,file=trim(filename),status='old')
write(filename,"('./Output/',A,'.dat')") trim(version)
open(101,file=filename)
write(*,*)
nLines=0;nDataPoints=0
do i=1,100 !Figure out how many lines and data points are in the file
  read(100,'(A)',end=50) tmp
  if(tmp(1:4).eq.':01:')then
    nDataPoints=nDataPoints+1
  else
    nLines=nLines+1
  end if
  ! if(i.eq.3)read(100,1001) date, time !Read the date and time of the flyby
  ! if(i.eq.3)write(*,1000) trim(version),date,time !Write to screen
  ! if(i.ge.16)read(100,1002) Jenergy(i-15),Jintensity(i-15)
end do
50 continue
!Allocate the variables based on the number of data points
allocate(Jenergy(nDataPoints),Jintensity(nDataPoints))
allocate(Jebins(nDataPoints+1),Jflux(nDataPoints),JebinWidth(nDataPoints))
nJebins=nDataPoints+1
Jenergy=0.0;Jintensity=0.0;Jebins=0.0;Jflux=0.0
rewind(100) !Rewind the file back to the beginning
!Read in the data points into variables
do i=1,100
  if(i.le.nLines)then
    read(100,*)
  else
    read(100,1002,end=51) Jenergy(i-nLines),Jintensity(i-nLines)
  end if
end do
51 continue
close(100) !Close JEDI measurement file
!Calculate the energy bin widths. This is variable, dependent upon how JEDI is
!collaborated. So rather than worrying about that, I assume the data points are
!in the middle of energy bins, and calculate them based on that.
Jebins(1)=JebinsMIN
Jebins(nDataPoints+1)=JebinsMAX
do i=2,nDataPoints
  Jebins(i)=Jenergy(i-1)+(Jenergy(i)-Jenergy(i-1))/2
  ! write(*,*) Jebins(i)-Jebins(i-1),Jenergy(i-1),Jebins(i),Jenergy(i)
end do
!Calculate the width of each energy bin
do i=1,nDataPoints
  JebinWidth(i)=Jebins(i+1)-Jebins(i)
end do
write(*,*)
write(*,*) trim(version)
if(spec.eq.oxy)then
  write(*,*) 'Oxygen output'
elseif(spec.eq.sul)then
  write(*,*) 'Sulfur output'
end if
write(*,1003)'Energy Bin:','Energy [keV/u]:','JEDI Intensity:',&
  'Energy Bin Width:','Normalized Flux:' !Write out general information
write(101,1003)'Energy Bin:','Energy [keV/u]:','JEDI Intensity:',&
  'Energy Bin Width:','Normalized Flux:' !Write out general information
do i=1,nDataPoints !Convert to [counts/cm^2/s]
!* The first 3 energy bins include both sulfur and oxygen. I'm assuming a 2:1
!* ratio of oxygen:sulfur (from SO_2)
  if(i.le.3)then
    if(spec.eq.oxy)then
      Jintensity(i)=Jintensity(i)*2/3
    elseif(spec.eq.sul)then
      Jintensity(i)=Jintensity(i)*1/3
    end if
  end if
  Jflux(i)=Jintensity(i)*2*pi*JebinWidth(i)
  write(*,1004)Jenergy(i),nint(Jenergy(i)/(16*spec)),Jintensity(i),&
               JebinWidth(i),Jflux(i)
  write(101,1004)Jenergy(i),nint(Jenergy(i)/(16*spec)),Jintensity(i),&
                 JebinWidth(i),Jflux(i)
end do
if(spec.eq.oxy)then
  write(101,*) 'Total number of oxygen ions:',sum(Jflux)
  write(101,*)
elseif(spec.eq.sul)then
  write(101,*) 'Total number of sulfur ions:',sum(Jflux)
  write(101,*)
end if
!**************************** Interpolate the data *****************************
!nInterp includes the JEDI energy data points as well
!nInterp also extrapolates two points lower and higher to cover the full range
!of the JEDI instrument.
nInterp=nJebins+(nJebins-1)*nSteps
allocate(Ienergy(nInterp),Iintensity(nInterp),Iebins(nInterp),IFlux(nInterp))
Ienergy=0.0;Iintensity=0.0;Iebins=0.0;Iflux=0.0
j=0
Ienergy(1)=log(JebinsMIN)
Ienergy(2)=Ienergy(1)+(log(Jenergy(1))-Ienergy(1))/2
Ienergy(nInterp-1)=log(JebinsMAX)-(log(JebinsMAX)-log(Jenergy(nDataPoints)))/2
Ienergy(nInterp)=log(JebinsMAX)
!Create interpolated energy points that include the original points (linearly)
do i=3,nInterp-2
  if(mod(i-3,nSteps+1).eq.0)j=j+1
  if(j.lt.nDataPoints)then
    Ienergy(i)=log(Jenergy(j))+& !Want a log interpolation
              (log(Jenergy(j+1))-log(Jenergy(j)))*&
              real(mod(i-3,nSteps+1))/real((nSteps+1))
  else
    Ienergy(i)=log(Jenergy(j))+& !Want a log interpolation
              (log(JebinsMAX)-log(Jenergy(j)))*&
              real(mod(i-3,nSteps+1))/real((nSteps+1))
  end if
end do
Ienergy=exp(Ienergy) !Put back into a normal number
!Loglog linearly interpolate the intensities
j=0
do i=3,nInterp-3
  if(mod(i-3,nSteps+1).eq.0)j=j+1
  if(Jintensity(j+1).lt.1e-8)Jintensity(j+1)=1.0e-20 !Make small rather than zero
  Iintensity(i)=log(Jintensity(j))+(log(Ienergy(i))-log(Jenergy(j)))*&
                (log(Jintensity(j+1))-log(Jintensity(j)))/&
                (log(Jenergy(j+1))-log(Jenergy(j)))
  Iintensity(i)=exp(Iintensity(i)) !Put back into a normal number
  if(Jintensity(j).lt.1.0e-5)Iintensity(i)=0.0
  if(Iintensity(i).lt.1.0e-5)Iintensity(i)=0.0
  ! if(isnan(Iintensity(i)))Iintensity(i)=0.0
  ! write(*,*) i,j,Iintensity(i),Jintensity(j),Jintensity(j+1)
end do
Iintensity(nInterp-2)=Jintensity(nDataPoints)
!Linearly extrapolate the first two data points
Iintensity(2)=Iintensity(4)+((Ienergy(2)-Ienergy(4))/(Ienergy(3)-Ienergy(4)))*&
              (Iintensity(3)-Iintensity(4))
Iintensity(1)=Iintensity(3)+((Ienergy(1)-Ienergy(3))/(Ienergy(2)-Ienergy(3)))*&
              (Iintensity(2)-Iintensity(3))
!Linearly extrapolate the final two data points
! Iintensity(nInterp-1)=Iintensity(nInterp-3)+&
!                       ((Ienergy(nInterp-1)-Ienergy(nInterp-3))/&
!                       (Ienergy(nInterp-2)-Ienergy(nInterp-3)))*&
!                       (Iintensity(nInterp-2)-Iintensity(nInterp-3))
Iintensity(nInterp-1)=log(Iintensity(nInterp-3))+&
                    ((log(Ienergy(nInterp-1))-log(Ienergy(nInterp-3)))/&
                    (log(Ienergy(nInterp-2))-log(Ienergy(nInterp-3))))*&
                    (log(Iintensity(nInterp-2))-log(Iintensity(nInterp-3)))
Iintensity(nInterp-1)=exp(Iintensity(nInterp-1))
Iintensity(nInterp)=log(Iintensity(nInterp-2))+&
                    ((log(Ienergy(nInterp))-log(Ienergy(nInterp-2)))/&
                    (log(Ienergy(nInterp-1))-log(Ienergy(nInterp-2))))*&
                    (log(Iintensity(nInterp-1))-log(Iintensity(nInterp-2)))
Iintensity(nInterp)=exp(Iintensity(nInterp))
!************************* Find new energy bin widths **************************
do i=1,nInterp !Assume the new energy points are the midpoint of the bins
  if(i.eq.1)then !Need to account for lowest energy value
    Iebins(i)=(Ienergy(i+1)-Ienergy(i))/2.0+(Ienergy(i)-JebinsMIN)
  elseif(i.eq.nInterp)then !Need to account for highest energy value
    Iebins(i)=(JebinsMAX-Ienergy(i))+(Ienergy(i)-Ienergy(i-1))/2.0
  else
    Iebins(i)=(Ienergy(i+1)-Ienergy(i))/2.0+(Ienergy(i)-Ienergy(i-1))/2.0
  end if
end do
write(*,*)
!************************ Calculate interpolated fluxes ************************
write(*,*)
if(spec.eq.oxy)then
  write(*,*) 'Interpolated Oxygen Values'
  write(101,*) 'Interpolated Oxygen Values'
elseif(spec.eq.sul)then
  write(*,*) 'Interpolated Sulfur Values'
  write(101,*) 'Interpolated Sulfur Values'
end if
write(*,1003)'Energy Bin:','Energy [keV/u]:','JEDI Intensity:',&
  'Energy Bin Width:','Normalized Flux:' !Write out general information
write(101,1003)'Energy Bin:','Energy [keV/u]:','JEDI Intensity:',&
  'Energy Bin Width:','Normalized Flux:' !Write out general information
do i=1,nInterp !Convert to [counts/cm^2/s]
  Iflux(i)=Iintensity(i)*2*pi*Iebins(i)
  write(*,1004)Ienergy(i),nint(Ienergy(i)/(16*spec)),Iintensity(i),&
               Iebins(i),Iflux(i)
  write(101,1004)Ienergy(i),nint(Ienergy(i)/(16*spec)),Iintensity(i),&
                 Iebins(i),Iflux(i)
end do
write(*,*) sum(Iebins*Iintensity),sum(JebinWidth*Jintensity)

if(spec.eq.oxy)then
  write(101,*) 'Total number of oxygen ions:',sum(Iflux)
elseif(spec.eq.sul)then
  write(101,*) 'Total number of sulfur ions:',sum(Iflux)
end if
close(101)

! 1000 format('FILE: ',A5,'.d2s',/,'DATE: ',A10,/,'TIME: ',A12)
! 1001 format(39X,A10,1X,A12)
1002 format(5X,ES12.9,1X,ES13.10)
1003 format(3x,A11,2x,A15,2x,A15,2x,A17,2x,A16)
1004 format(F11.3,1x,I12,5x,F15.4,2x,F17.2,2x,F16.3)

end subroutine
! end program
