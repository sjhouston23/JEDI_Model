subroutine XRayOpacity(NoOpacSpec,WaveL,nLines,vAngle,OpacSpec)
!*******************************************************************************
!* Created by Stephen J. Houston 1.25.18
!*******************************************************************************
!* This subroutine takes the synthetic ion X-ray spectra produced by
!* XRaySpectra.f08 and the absorption cross-sections interpolated by
!* Absorption.f08 and generates X-ray spectra that include opacity effects.
!*******************************************************************************
!*
!* 	Input:
!*		NoOpacSpec --> Initial X-ray spectrum with no opacity effects
!*		 Type: 2-dimensional, real array (intensity vs. wavelength and altitude)
!*		 Units: photons/cm^3/s
!*
!*		WaveL --> Wavelength of emission (in eV) and photon intensity
!*		 Type: 2-dimensional, real array (1-Photon energy, 2-Intensity)
!*		 Units: eV, normalized photons/ion
!*
!*		nLines --> Number of emission lines
!*		 Type: Integer
!*		 Units: None
!*
!*    vAngle --> Viewing angle for opacity effets
!*		 Type: Integer
!*		 Units: Degrees
!*
!*   Returns:
!*    OpacSpec --> X-Ray spectrum with opacity effects
!*		 Type: 2-dimensional, real array (Wavelength for each altitude)
!*		 Units: photons/cm^3/s
!*
!*******************************************************************************

implicit real*8(a-h,o-z)

!*******************************************************************************
integer AtmosLen
parameter(AtmosLen=1544)

real*8,intent(in) :: NoOpacSpec(AtmosLen,nLines),WaveL(2,nLines),vAngle
real*8,intent(out) :: OpacSpec(AtmosLen,nLines)

integer H2,He,CH4

parameter(Rj=71492.0e5) !Radius of Jupiter [km]
parameter(pi=4.0d0*atan(1.0d0))
parameter(h=4.135668e-15) !Planck constant [ev s]
parameter(c=2.9979e17) !Speed of light [nm/s]
parameter(H2=1,He=2,CH4=3)

integer Alt,WL !Altitude, wavelength

real*8 AbsXS(3)
real*8,dimension(AtmosLen) :: Chapman
real*8,dimension(AtmosLen),save :: Altitude,cdH2,cdHe,cdCH4,cdH,cdTot,SH,Opacity

!*************************** Read in the Atmosphere ****************************
if(Altitude(1).gt.1) goto 10
open(unit=100,file='./Atmosphere/JupiterColumnDensity_2km.dat',status='old')
read(100,*) !Skip header lines
do Alt=1,AtmosLen !Read in the column density values
  read(100,*) Altitude(Alt),cdH2(Alt),cdHe(Alt),cdCH4(Alt),cdH(Alt),cdTot(Alt),&
              SH(Alt) !SH - Scale height in cm
  !Altitude(1) is the top of the atmosphere, 2998.0 km; Altitude(1544)=-88.0 km
end do
close(100) !Close column density file
10 continue
!******************************** Main Program *********************************
if(vAngle.lt.81.0)then !If the viewing angle is below 81° then we want sec()
  do Alt=1,AtmosLen
    Chapman(Alt)=1/cos(vAngle*pi/180.0) !Chapman function for <80°
  end do
else !If it's greater than 81°, then we want to use the chapman function
  Chapman=sqrt((Rj/SH)*pi/2) !Chapman function for 90°
end if
do WL=1,nLines !Loop through each line emission wavelength
  WaveLnm=h*c/WaveL(1,WL) !Convert eV to nm
  call Absorption(WaveLnm,AbsXS) !Get the absorption cross-sections vs. wl
  do Alt=1,AtmosLen !Loop through the entire atmosphere
    TauH2=2*cdH2(Alt)*AbsXS(H2)*Chapman(Alt) !AbsXS is for H
    TauHe=cdHe(Alt)*AbsXS(He)*Chapman(Alt)
    TauCH4=cdCH4(Alt)*AbsXS(CH4)*Chapman(Alt) !AbsXS is for C
    Opacity(Alt)=exp(-(TauH2+TauHe+TauCH4))
    OpacSpec(Alt,WL)=NoOpacSpec(Alt,WL)*Opacity(Alt)
  end do
end do

1004 format(1x,F8.2,50(1x,ES9.3E2))

end subroutine
