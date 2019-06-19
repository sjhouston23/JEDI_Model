subroutine combCXDE(version,vAngle)
!*******************************************************************************
!* Created by Stephen J. Houston 3.5.18
!*******************************************************************************
!* This program combines the charge exchange and direct excitation emission
!* lines (points) for both sulfur and oxygen.
!*******************************************************************************

implicit real*8(a-h,o-z)

!*******************************************************************************
integer ChS,Ang
parameter(MaxnLines=1000,nChS=9,nProc=2,MaxWL=500)
integer nLines(nProc,nChS)

real*8 vAngle(3)
real*8,dimension(nProc,nChS,MaxWL) :: WL,Flux
real*8,dimension(nChS,MaxWL) :: TotalFlux

character(len=100) :: filename,version
!******************************** Main Program *********************************
do Ang=1,4
  !*** Initialize Variables
  Flux=0.0;WL=0.0;nLines=0.0;TotalFlux=0.0
  do ChS=1,nChS !O6+, O7+, S8+ - S14+
    if(Ang.le.3)then
      if(ChS.le.2)write(filename,"('./Output/XRay/',A,'/Points/',I0,'/Oxy_',I0,&
        '_CX.dat')") trim(version),nint(vAngle(Ang)),ChS+5
      if(ChS.ge.3)write(filename,"('./Output/XRay/',A,'/Points/',I0,'/Sul_',I0,&
        '_CX.dat')") trim(version),nint(vAngle(Ang)),ChS+5
    else if(Ang.eq.4)then
      if(ChS.le.2)write(filename,"('./Output/XRay/',A,'/Points/Oxy_',I0,&
        '_CX.dat')") trim(version),ChS+5
      if(ChS.ge.3)write(filename,"('./Output/XRay/',A,'/Points/Sul_',I0,&
        '_CX.dat')") trim(version),ChS+5
    end if
    open(100,file=filename,status='old')
    do i=1,MaxnLines
      read(100,*,end=500) WL(1,ChS,i),Flux(1,ChS,i) !Read in the file
      nLines(1,ChS)=nLines(1,ChS)+1
    end do
    500 continue
    close(100)
    if(Ang.le.3)then
      if(ChS.le.2)write(filename,"('./Output/XRay/',A,'/Points/',I0,'/Oxy_',I0,&
        '_DE.dat')") trim(version),nint(vAngle(Ang)),ChS+5
      if(ChS.ge.3)write(filename,"('./Output/XRay/',A,'/Points/',I0,'/Sul_',I0,&
        '_DE.dat')") trim(version),nint(vAngle(Ang)),ChS+5
    else if(Ang.eq.4)then
      if(ChS.le.2)write(filename,"('./Output/XRay/',A,'/Points/Oxy_',I0,&
        '_DE.dat')") trim(version),ChS+5
      if(ChS.ge.3)write(filename,"('./Output/XRay/',A,'/Points/Sul_',I0,&
        '_DE.dat')") trim(version),ChS+5
    end if
    open(101,file=filename,status='old')
    do i=1,MaxnLines
      read(101,*,end=501) WL(2,ChS,i),Flux(2,ChS,i) !Read in the file
      nLines(2,ChS)=nLines(2,ChS)+1
    end do
    501 continue
    close(101)
    j=1
    do i=1,nLines(1,ChS)
      if(WL(1,ChS,i).eq.WL(2,ChS,j))then
        TotalFlux(ChS,i)=Flux(1,ChS,i)+Flux(2,ChS,j)
        j=j+1
      else
        TotalFlux(ChS,i)=Flux(1,ChS,i)
      end if
      if(i.eq.nLines(1,ChS).and.j-1.ne.nLines(2,ChS))&
        write(*,*)'Line mismatch for: ',ChS+5,i,nLines(1,ChS),j-1,nLines(2,ChS)
    end do
  end do

  do ChS=1,nChS
    if(Ang.le.3)then
      if(ChS.le.2)write(filename,"('./Output/XRay/',A,'/Points/',I0,&
        '/CX+DE_Oxy_',I0,'.dat')") trim(version),nint(vAngle(Ang)),ChS+5
      if(ChS.ge.3)write(filename,"('./Output/XRay/',A,'/Points/',I0,&
        '/CX+DE_Sul_',I0,'.dat')") trim(version),nint(vAngle(Ang)),ChS+5
    else if(Ang.eq.4)then
      if(ChS.le.2)write(filename,"('./Output/XRay/',A,'/Points/CX+DE_Oxy_',I0,&
        '.dat')") trim(version),ChS+5
      if(ChS.ge.3)write(filename,"('./Output/XRay/',A,'/Points/CX+DE_Sul_',I0,&
        '.dat')") trim(version),ChS+5
    end if
    open(102,file=filename)
    do i=1,nLines(1,ChS)
      write(102,1000) WL(1,ChS,i),TotalFlux(ChS,i)
    end do
  end do
end do

1000 format(1x,F8.3,2x,ES9.3E2)

end subroutine
