program JEDIModel
!*******************************************************************************
!* Created by Stephen J. Houston 5.29.19
!*******************************************************************************
!* This program calculates all of the results (X-ray, secondary electron
!* production, ionization rates, etc...) given an input JEDI ion flux, both
!* oxygen and sulfur.
!*******************************************************************************

implicit real*8(a-h,o-z)

!**************************** Variable Declaration *****************************
real*8 Angles(3)
character(len=1) ans
character(len=100) version,check
logical ex
!****************************** Data Declaration *******************************
data Angles/0.0,60.0,80.0/ !Can only do 3 at a time, but this could be changed.
!*******************************************************************************
write(*,*)
write(*,*) "What is the name of the JEDI ion spectrum file you wish to open?"
write(*,*) "Don't include the extension in the file name (e.g. if you want to &
 analyze PJ7-1.ds2, type PJ7-1)."
read(*,*) version
!Oxygen JEDI calculations
!First check to see if the files already exists
write(check,"('./Output/Oxy/',A,'/Total_Photon_Prod.dat')") trim(version)
inquire(file=check,exist=ex)
if(ex)then
  write(*,*)
  write(*,*) 'Output oxygen files already exists for this JEDI pass.'
  write(*,*) 'Using these files will not allow for the calculation of the S:O &
             &ratio.'
  write(*,*) 'Do you want to use the existing files? (Y/N)'
  read(*,*) ans
  if(ans.eq.'Y')then
    write(*,*)
    write(*,*) 'Skipping to sulfur precipiation...'
    goto 100
  end if
  write(*,*)
  write(*,*) 'Recalculating the interpolations and generating new output &
             &files...'
  call sleep(2)
end if
call JEDIOxyPrecip(version,Oxy) !Returns total number of oxygen ions
100 continue
!Sulfur JEDI calculations
!First check to see if the files already exists
write(check,"('./Output/Sul/',A,'/Photons_Total.dat')") trim(version)
inquire(file=check,exist=ex)
if(ex)then
  write(*,*)
  write(*,*) 'Output sulfur files already exists for this JEDI pass.'
  write(*,*) 'Using these files will not allow for the calculation of the S:O &
             &ratio.'
  write(*,*) 'Do you want to use the existing files? (Y/N)'
  read(*,*) ans
  if(ans.eq.'Y')then
    write(*,*)
    write(*,*) 'Skipping to X-ray calculations...'
    goto 101
  end if
  write(*,*)
  write(*,*) 'Recalculating the interpolations and generating new output &
             &files...'
  call sleep(2)
end if
call JEDISulPrecip(version,Sul) !Returns total number of sulfur ions
write(*,*)
write(*,1000) 'This spectrum has an S:O ratio of -',Sul/Oxy,'(',Sul,':',Oxy,')'
101 continue
!X-Ray calculations
write(*,*)
write(*,*) 'Calculating X-ray production from charge exchange (CX) processes &
           &for the original atmosphere...'
call XRaySpectraCX(version,Angles)
write(*,*)
write(*,*) 'Calculating X-ray production from direct excitation (DE) processes &
           &for the original atmosphere...'
call XRaySpectraDE(version,Angles)
write(*,*)
write(*,*) 'Calculating X-ray production from charge exchange (CX) processes &
           &for the well-mixed atmosphere...'
call XRaySpectraCXAtm2(version,Angles)
write(*,*)
write(*,*) 'Calculating X-ray production from direct excitation (DE) processes &
           &for the well-mixed atmosphere...'
call XRaySpectraDEAtm2(version,Angles)
!Combine the X-Ray line emission points from CX and DE.
!This is useful when sending data to Will Dunn so he can plug it into the
!instrument response functions for comparing spectra.
call combCXDE(version,Angles)
call combCXDEAtm2(version,Angles)
write(*,*)
write(*,*) 'DONE'
write(*,*)

1000 format(A36,F5.2,2x,A1,F11.2,1x,A1,F11.2,A1)

end program
