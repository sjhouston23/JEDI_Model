subroutine Absorption(WaveLI,AbsXSI)
!*******************************************************************************
!* Created by Stephen J. Houston 1.22.18
!*******************************************************************************
!* This subroutine calculates the absorption cross-section for H, He, and C
!* (from CH4) for the calculation of the opacity effects on X-ray production
!*******************************************************************************
!*
!* 	Input:
!*		WaveLI --> Wavelength of the photon
!*		 Type: Real
!*		 Units: nm
!*
!*   Returns:
!*    AbsXSI --> Interpolated value of the absorption cross-section
!*		 Type: 3 component, 1-D real array
!*		 Units: cm^2s
!*
!*******************************************************************************

implicit real*8(a-h,o-z)

!**************************** Variable Declaration *****************************
real*8,intent(in) :: WaveLI
real*8,intent(out) :: AbsXSI(3)

integer Spec,H,He,C

parameter(nMaxLines=100) !Max number of lines for absorption cross-section files
parameter(nSpec=3) !Number of species; H, He, C (from CH4)
parameter(H=1,He=2,C=3)

integer,save :: nLines(nSpec)

real*8 AbsXSI_tmp
real*8,allocatable,dimension(:) :: WaveLVec,AbsXSVec,AbsXSVec2
real*8,allocatable,dimension(:,:),save :: WaveL,ScatXS,AbsXS

!******************************** Read in XS's *********************************
if(nLines(1).gt.0) goto 30
open(unit=100,file='./AbsorptionXS/H_AbsXS.dat',status='old') !Hydrogen XS
open(unit=101,file='./AbsorptionXS/He_AbsXS.dat',status='old') !Helium XS
open(unit=102,file='./AbsorptionXS/C_AbsXS.dat',status='old') !Carbon (CH4) XS
read(100,*);read(101,*);read(102,*)
nLines=0
do Spec=1,nSpec !Loop through each species
  do j=1,nMaxLines !Loop through to see how many lines each file has
    read(99+Spec,*,end=10) !Get out once the end is reached
    nLines(Spec)=nLines(Spec)+1
  end do !End max lines do-loop
10 continue
  rewind(99+Spec) !Rewind the files
end do !End each species do-loop
allocate(WaveL(nSpec,maxval(nLines))) !Wavelength [nm]
allocate(ScatXS(nSpec,maxval(nLines))) !Scattering cross-section [cm^2]
allocate(AbsXS(nSpec,maxval(nLines))) !Absorption cross-section [cm^2]
WaveL=0.0;ScatXS=0.0;AbsXS=0.0 !Initialize all to zero
read(100,*);read(101,*);read(102,*) !Skip the header line
do Spec=1,nSpec !Loop through each species and file
  do j=nLines(Spec),1,-1!1,nLines(Spec) !Loop through every line
    read(99+Spec,*) WaveL(Spec,j),ScatXS(Spec,j),AbsXS(Spec,j) !Read in the XS's
  end do !End nlines do-loop
  close(99+Spec)
end do !End each species/file do-loop
30 continue
!******************************** Main Program *********************************
!*** Interpolate the cross-sections for a specific wavelength:
do Spec=1,nSpec
  allocate(WaveLVec(nLines(Spec)),AbsXSVec(nLines(Spec)),&
           AbsXSVec2(nLines(Spec)))
  WaveLVec=0.0;AbsXSVec=0.0;AbsXSVec2=0.0 !Initialize all to zero
  do j=1,nLines(Spec)
    WaveLVec(j)=WaveL(Spec,j) !Create a wavelength vector for spline
    AbsXSVec(j)=AbsXS(Spec,j) !Create an absorption XS vector for spline
  end do
  call spline(log(WaveLVec),log(AbsXSVec),nLines(Spec),AbsXSVec2)
  k=1
20 continue
  if(Spec.eq.C.and.WaveLI.ge.WaveLVec(k+1))then
    k=k+1
    goto 20
  end if
  call splineinterp(log(WaveLI),log(WaveLVec),log(AbsXSVec),nLines(Spec),&
                    AbsXSVec2,AbsXSI_tmp)
  if(Spec.eq.C.and.WaveLI.ge.1.and.WaveLI.le.10)then !linearly interpolate
    AbsXSI_tmp=log(AbsXSVec(k))+(log(WaveLI)-log(WaveLVec(k)))*&
    (log(AbsXSVec(k+1))-log(AbsXSVec(k)))/&
    (log(WaveLVec(k+1))-log(WaveLVec(k)))
  end if
  AbsXSI(Spec)=exp(AbsXSI_tmp)
  deallocate(WaveLVec,AbsXSVec,AbsXSVec2)
end do

end subroutine
