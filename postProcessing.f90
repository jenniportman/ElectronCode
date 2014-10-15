program postprocessing
implicit none

character(80) :: fileIn,fileOut, fileOut2, fileOut3, fileOut4, fileOut5
integer :: nStep, nProc, nMacroParticle, nFreq
integer :: iMode

write(*,*) "What do you want to do?"
write(*,*) "1=AG parameters, 2=density plot, 3=emittance vs time, 4=radial profile, 5=longitudinal profile"
read(*,*) iMode

write(*,*) "Insert total number of processors"
read(*,*) nProc
write(*,*) "Insert number of electrons in each macroparticle"
read(*,*) nMacroParticle

if(iMode == 2 .or. iMode == 4 .or. iMode == 5) then
    write(*,*) "Insert time step number"
    read(*,*) nStep
else
    write(*,*) "Insert total number of time steps"
    read(*,*) nStep
    write(*,*) "Insert file save frequency (1 or 10)"
    read(*,*) nFreq
end if

if(iMode == 1) then
    write(fileOut,*) "AGParameters.txt"
    open(unit=41,file=trim(fileOut),status="unknown",recl=600)
else if(iMode == 2) then
    write(fileOut,'(i,a)') nStep,"-xpx.txt"
    open(unit=41,file=trim(fileOut),status="unknown", recl=500)
    write(fileOut2,'(i,a)') nStep,"-zpz.txt"
    open(unit=42,file=trim(fileOut2),status="unknown", recl=500)
    write(fileOut3,'(i,a)') nStep,"-vector.txt"
    open(unit=43,file=trim(fileOut3),status="unknown", recl=500)
    write(fileOut4,'(i,a)') nStep,"-xz.txt"
    open(unit=44,file=trim(fileOut4),status="unknown", recl=500)
    write(fileOut5,'(i,a)') nStep,"-velocity.txt"
    open(unit=45,file=trim(fileOut5),status="unknown", recl=500)
else if(iMode == 3) then
    write(fileOut,*) "emittance_data_norm.txt"
    open(unit=41,file=trim(fileOut),status="unknown",recl=400)
else if(iMode == 4) then
    write(fileOut,'(i,a)') nStep,"-radial_profile.txt"
    open(unit=41,file=trim(fileOut),status="unknown",recl=400)
else if(iMode == 5) then
    write(fileOut,'(i,a)') nStep,"-long_profile.txt"
    open(unit=41,file=trim(fileOut),status="unknown",recl=400)
end if

write(fileIn,'(a)') "screen.txt"
open(unit=50,file=trim(fileIn),status="old")

if(iMode == 1) then
    call AGParameters(nStep,nFreq,nMacroParticle,nProc)
else if(iMode == 2) then
    call densityPlot(nStep,nMacroParticle,nProc)
else if(iMode == 3) then
    call emittance(nStep,nFreq,nMacroParticle,nProc)
else if(iMode == 4) then
    call radialProfile(nStep,nProc)
else if(iMode == 5) then
    call longitudinalProfile(nStep,nProc)
end if

close(41)
close(50)

if(iMode == 2) then
    close(42)
    close(43)
    close(44)
    close(45)
end if

contains

!-------------------------------------------------------------------
! subroutine AGParameters
! Calculates parameters for Analytic Gaussian model from photoemission simulation
!-------------------------------------------------------------------
subroutine AGParameters(nStep,nFreq,nMacroParticle,nProc)
implicit none
integer, intent(in) :: nStep, nFreq, nMacroParticle, nProc

integer :: i,stat,j,k
real, dimension(3) :: meanX,meanX2,meanP,meanP2,meanXP,xVect,pVect
real, dimension(3) :: sigmaX, gammaXP, etaP, etaP0
real, allocatable :: time(:)
real :: x1, x2, vz, gammaL
real :: position_conv=1.0e6, time_conv=1e12, momentum_conv, speed_light
character(80) :: fileIn

speed_light = 299792458*1e-6
momentum_conv = speed_light/(0.51099892811*nMacroParticle)

write(41,*) "time (ps), Ne, z (um), vz (um/ps), sigma_x,y,z^1/2 (um), gamma_x,y,z (um*(me*um/ps)), stand dev p ((me*um/ps))"

allocate(time(nStep))
do i=1,nStep
    read(50,*) x1, x2, time(i)
end do

time(:) = time(:)*time_conv
do i=1*nFreq,nStep, nFreq
    meanX=0.0
    meanX2=0.0
    meanP=0.0
    meanP2=0.0
    meanXP=0.0
    j=0

    do k=1,nProc
        stat=0
        ! Generate file name and open it
        if (k < 10) then
            write(fileIn,'(i,a,i1,a)') i,"-x-",k,".txt"
        else
            write(fileIn,'(i,a,i2,a)') i,"-x-",k,".txt"
        end if
        open(unit=40,file=trim(fileIn),status="old")

        do while(stat==0)
            j=j+1
            read(40,*,iostat=stat) xVect,pVect
            xVect(:) = xVect(:)*position_conv
            pVect(:) = pVect(:)*momentum_conv
            if(stat == 0) then
                meanX(:)=meanX(:)+xVect(:)
                meanX2(:)=meanX2(:)+xVect(:)**2.0
                meanP(:)=meanP(:)+pVect(:)
                meanP2(:)=meanP2(:)+pVect(:)**2.0
                meanXP(:)=meanXP(:)+xVect(:)*pVect(:)
            end if
        end do
        close(40)
        j=j-1
    end do

    write(*,*) "Nstep", i, "Nel", j
    sigmaX(:)=(real(j)*meanX2(:)-meanX(:)*meanX(:))/(real(j)*real(j))
    gammaXP(:)=(real(j)*meanXP(:)-meanX(:)*meanP(:))/(real(j)*real(j))
    etaP(:)=(real(j)*meanP2(:)-meanP(:)*meanP(:))/(real(j)*real(j))
    etaP0(:) = etaP(:)-gammaXP(:)**2.0/sigmaX(:)
    meanX(:)=meanX(:)/real(j)
    meanP(:)=meanP(:)/real(j)
    gammaL=sqrt(1+(meanP(3)/speed_light)**2.0)
    vz=meanP(3)/gammaL

    write(41,*) time(i),meanX(3),vz,gammaL,sqrt(sigmaX(1)),sqrt(sigmaX(3)),gammaXP(1),gammaXP(3),sqrt(etaP0(1)),sqrt(etaP0(3)),j*nMacroParticle
    !write(41,*) time(i),j*nMacroParticle,meanX(3),vz,sqrt(sigmaX),gammaXP,sqrt(etaP)
end do
write(*,*) "initial values to use"
write(*,*) "Ne", j*nMacroParticle
write(*,*) "velocity vz",vz
write(*,*) sqrt(sigmaX(1)),sqrt(sigmaX(3))*gammaL,gammaXP(1),gammaXP(3),sqrt(etaP0(1)),sqrt(etaP0(3))/gammaL
!write(*,*) sqrt(sigmaX),gammaXP,sqrt(etaP0)
write(*,*) "length acceleration gap", 45.0E3-meanX(3)

end subroutine AGParameters

!-------------------------------------------------------------------
! subroutine densityPlot
! Calculates data for density plot
!-------------------------------------------------------------------
subroutine densityPlot(nStep,nMacroParticle,nProc)
implicit none
integer, intent(in) :: nStep, nMacroParticle, nProc

real, parameter :: nXSigma=3.0, nZSigma=4.0, nPzSigma=3.0, nPxSigma=3.0
integer, parameter :: nGrid=160
integer :: i,stat,j,k
integer :: ix, iz, ipz, ipx 
integer, dimension(nGrid+1, nGrid+1) :: densityArrayXZ, densityArrayZPZ, densityArrayXPX, nParticleArray
real, dimension(nGrid+1, nGrid+1) :: vectorDx, vectorDz
real, dimension(nGrid+1, nGrid+1) :: densityArrayXZSmooth,densityArrayZPZSmooth, densityArrayXPXSmooth
real, dimension(nGrid+1, nGrid+1) :: velocityXArray, velocity2XArray, stdVelocityXArray, velocityZArray, velocity2ZArray, stdVelocityZArray
real, dimension(15) :: inputVector
real, dimension(6) :: particleVector
real :: zMin, zMax, xMin, xMax, pzMin, pzMax, pxMin, pxMax, dzGrid, dxGrid, dpzGrid, dpxGrid
real :: meanPx, sigmaPx, meanPz, sigmaPz, meanX, meanZ
real :: pxMaxV, pzMaxV, stdXMax, stdZMax
character(80) :: fileIn

do i=1,nStep-1
    read(50,*)
end do
read(50,*) inputVector

meanX=inputVector(4)
meanZ=inputVector(6)
meanPx=inputVector(10)/nMacroParticle
meanPz=inputVector(12)/nMacroParticle
sigmaPx=inputVector(13)/nMacroParticle
sigmaPz=inputVector(15)/nMacroParticle

! Calculate grid ranges and spacing
!zMin=0.0
zMin=meanZ-nZSigma*inputVector(9)
zMax=meanZ+nZSigma*inputVector(9)
dzGrid=(zMax-zMin)/real(nGrid)

xMin=meanX-nXSigma*inputVector(7)
xMax=meanX+nXSigma*inputVector(7)
dxGrid=(xMax-xMin)/real(nGrid)

pzMin=meanPz-nPzSigma*sigmaPz
pzMax=meanPz+nPzSigma*sigmaPz
dpzGrid=(pzMax-pzMin)/real(nGrid)

pxMin=meanPx-nPxSigma*sigmaPx
pxMax=meanPx+nPxSigma*sigmaPx
dpxGrid=(pxMax-pxMin)/real(nGrid)

write(*,*) "grid ranges"
write(*,*) "z from", zMin, "to", zMax
write(*,*) "x from", xMin, "to", xMax
write(*,*) "pz from", pzMin, "to", pzMax
write(*,*) "px from", pxMin, "to", pxMax
write(*,*) "nGrid",nGrid

densityArrayXZ=0
densityArrayZPZ=0
densityArrayXPX=0
velocityXArray=0.0
velocity2XArray=0.0
velocityZArray=0.0
velocity2ZArray=0.0
nParticleArray=0
vectorDx=0
vectorDz=0
k=0

do i=1,nProc
    ! Generate file name and open it
	if (i < 10) then
		write(fileIn,'(i,a,i1,a)') nStep,"-x-",i,".txt"
	else
		write(fileIn,'(i,a,i2,a)') nStep,"-x-",i,".txt"
	end if
	write(*,*) "reading file",trim(fileIn)
	open(unit=40,file=trim(fileIn),status="old")

    ! Start reading file and bin data
	stat=0
	do while(stat==0)
		k=k+1
		read(40,*,iostat=stat) particleVector
        particleVector(4:6) = particleVector(4:6)/nMacroParticle
        if(stat == 0) then
            ! density plot stuff
            ! z,x + vectors
            if(particleVector(3) < zMax .and. particleVector(3) > zMin) then
                iz=int((particleVector(3)-zMin)/dzGrid)+1
                if(particleVector(1) < xMax .and. particleVector(1) > xMin) then
                    ix=int((particleVector(1)-xMin)/dxGrid)+1
                    densityArrayXZ(iz,ix)=densityArrayXZ(iz,ix)+1
                    vectorDx(iz,ix)=vectorDx(iz,ix)+particleVector(4)
                    vectorDz(iz,ix)=vectorDz(iz,ix)+particleVector(6)
                    velocityXArray(iz,ix)=velocityXArray(iz,ix)+particleVector(4)
                    velocityZArray(iz,ix)=velocityZArray(iz,ix)+(particleVector(6)-meanPz)
                    velocity2XArray(iz,ix)=velocity2XArray(iz,ix)+particleVector(4)**2
                    velocity2ZArray(iz,ix)=velocity2ZArray(iz,ix)+(particleVector(6)-meanPz)**2
                    !velocityArray(iz,ix)=velocityArray(iz,ix)+sqrt(particleVector(4)**2+particleVector(6)**2)
                    !velocitySqArray(iz,ix)=velocitySqArray(iz,ix)+(particleVector(4)**2+particleVector(6)**2)
                    nParticleArray(iz,ix)=nParticleArray(iz,ix)+1
                end if
            ! z,pz
                if(particleVector(6) < pzMax .and. particleVector(6) > pzMin) then
                    ipz=int((particleVector(6)-pzMin)/dpzGrid)+1
                    densityArrayZPZ(iz,ipz)=densityArrayZPZ(iz,ipz)+1
                end if
            end if
            ! x,px
            if(particleVector(1) < xMax .and. particleVector(1) > xMin) then
                ix=int((particleVector(1)-xMin)/dxGrid)+1
                if(particleVector(4) < pxMax .and. particleVector(4) > pxMin) then
                    ipx=int((particleVector(4)-pxMin)/dpxGrid)+1
                    densityArrayXPX(ix,ipx)=densityArrayXPX(ix,ipx)+1
                end if
            end if 
		end if
	end do
	k=k-1
	close(40)
end do

call smoothing(densityArrayXZ,densityArrayXZSmooth) 
call smoothing(densityArrayZPZ,densityArrayZPZSmooth) 
call smoothing(densityArrayXPX,densityArrayXPXSmooth) 

do i=1,nGrid+1
    do j=1,nGrid+1
        if(nParticleArray(i,j) /= 0) then
            stdVelocityXArray(i,j) = sqrt(velocity2XArray(i,j)*real(nParticleArray(i,j))-(velocityXArray(i,j))**2.0)/real(nParticleArray(i,j))
            stdVelocityZArray(i,j) = sqrt(velocity2Zarray(i,j)*real(nParticleArray(i,j))-(velocityZArray(i,j))**2.0)/real(nParticleArray(i,j))
            !write(*,*) "TEST", vectorDz(i,j)**2,velocityzarray(i,j)
            !stdVelocityArray(i,j) = sqrt(velocityArray(i,j)*real(nParticleArray(i,j))-(vectorDz(i,j)+vectorDx(i,j))**2.0)/real(nParticleArray(i,j))
        end if
    end do
end do
! Calculate density arrays and write to file
pxMaxV=maxval(vectorDx)
pzMaxV=maxval(vectorDz)
stdXMax=maxval(stdVelocityXArray)
stdZMax=maxval(stdVelocityZArray)
do i=1,nGrid+1
    do j=1,nGrid+1
        write(41,*) xMin+(i-1)*dxGrid-meanX, pxMin+(j-1)*dpxGrid-meanPX, densityArrayXPXSmooth(i,j)
        write(42,*) zMin+(i-1)*dzGrid-meanZ, pzMin+(j-1)*dpzGrid-meanPZ, densityArrayZPZSmooth(i,j)
        write(44,*) zMin+(i-1)*dzGrid-meanZ, xMin+(j-1)*dxGrid-meanX, densityArrayXZSmooth(i,j)
        if(vectorDz(i,j) /= 0 .or. vectorDx(i,j) /= 0) then 
            write(43,*) zMin+(i-1)*dzGrid-meanZ, xMin+(j-1)*dxGrid-meanX, (vectorDz(i,j)-meanPz)/pzMaxV*dzGrid, vectorDx(i,j)/pxMaxV*dxGrid
        else
            write(43,*) zMin+(i-1)*dzGrid-meanZ, xMin+(j-1)*dxGrid-meanX, (vectorDz(i,j))/pzMaxV*dzGrid, vectorDx(i,j)/pxMaxV*dxGrid
        end if
        write(45,*) zMin+(i-1)*dzGrid-meanZ, xMin+(j-1)*dxGrid-meanX, stdVelocityXArray(i,j)/stdXMax, stdVelocityZArray(i,j)/stdZMax
    end do
    write(41,*)
    write(42,*)
    write(43,*)
    write(44,*)
    write(45,*)
end do

write(*,*) "rows read",k

end subroutine densityPlot


subroutine smoothing(ArrayIn,ArrayOut)
implicit none
integer, intent(in) :: ArrayIn(:,:)
real, intent(out) :: ArrayOut(:,:)

integer :: i, j, nGrid, densityCount
real :: densitySum, maxValue

nGrid = size(ArrayIn,1)

! Smooth arrays for nicer plots
do i=1,nGrid+1
    do j=1,nGrid+1
        densitySum=2*ArrayIn(i,j)
        densityCount=2
        if(i > 1) then
            densitySum=densitySum+ArrayIn(i-1,j)
            densityCount=densityCount+1
            if(j > 1) then
                densitySum=densitySum+ArrayIn(i-1,j-1)
                densityCount=densityCount+1
            end if
            if(j < nGrid+1) then
                densitySum=densitySum+ArrayIn(i-1,j+1)
                densityCount=densityCount+1
            end if
        end if
        if(i < nGrid+1) then
            densitySum=densitySum+ArrayIn(i+1,j)
            densityCount=densityCount+1
            if(j > 1) then
                densitySum=densitySum+ArrayIn(i+1,j-1)
                densityCount=densityCount+1
            end if
            if(j < nGrid+1) then
                densitySum=densitySum+ArrayIn(i+1,j+1)
                densityCount=densityCount+1
            end if
        end if
        if(j > 1) then
            densitySum=densitySum+ArrayIn(i,j-1)
            densityCount=densityCount+1
        end if
        if(j < nGrid+1) then
            densitySum=densitySum+ArrayIn(i,j+1)
            densityCount=densityCount+1
        end if
        ArrayOut(i,j)=real(densitySum)/real(densityCount)
    end do
end do
maxValue=maxval(ArrayOut)
ArrayOut = ArrayOut/maxValue
end subroutine smoothing

!-------------------------------------------------------------------
! subroutine emittance
! Calculates emittance as a function of time
!-------------------------------------------------------------------
subroutine emittance(nStep,nFreq,nMacroParticle,nProc)
implicit none
integer, intent(in) :: nStep, nFreq, nMacroParticle, nProc

character(80) :: fileIn
integer :: i,j,k,stat 
real, dimension(3) :: meanX,meanX2,meanP,meanP2,meanXP,rmsEmittance,xVect,pVect
real, dimension(3) :: sigmaX,sigmaP,sigmaXP
real, allocatable :: time(:)
real :: x1, x2

allocate(time(nStep))
do i=1,nStep
    read(50,*) x1, x2, time(i)
end do

do i=nFreq,nStep,nFreq
    meanX=0.0
    meanX2=0.0
    meanP=0.0
    meanP2=0.0
    meanXP=0.0
	j=0

    do k=1,nProc
        stat=0
        ! Generate file name and open it
        if (k < 10) then
            write(fileIn,'(i,a,i1,a)') i,"-x-",k,".txt"
        else
            write(fileIn,'(i,a,i2,a)') i,"-x-",k,".txt"
        end if
        open(unit=40,file=trim(fileIn),status="old")

        do while(stat==0)
            j=j+1
            read(40,*,iostat=stat) xVect,pVect
            if(stat == 0) then
                meanX(:)=meanX(:)+xVect(:)
                meanX2(:)=meanX2(:)+xVect(:)**2.0
                meanP(:)=meanP(:)+pVect(:)
                meanP2(:)=meanP2(:)+pVect(:)**2.0
                meanXP(:)=meanXP(:)+xVect(:)*pVect(:)
            end if
        end do
        close(40)
        j=j-1
    end do

    write(*,*) "Nel", j
    meanP(:)=meanP(:)/nMacroParticle
    meanP2(:)=meanP2(:)/nMacroParticle**2.d0
    meanXP(:)=meanXP(:)/nMacroParticle
    sigmaX(:)=(real(j)*meanX2(:)-meanX(:)*meanX(:))/(real(j)*real(j))
    sigmaP(:)=(real(j)*meanP2(:)-meanP(:)*meanP(:))/(real(j)*real(j))
    sigmaXP(:)=(real(j)*meanXP(:)-meanX(:)*meanP(:))/(real(j)*real(j))
	rmsEmittance(:)=sqrt(sigmaX(:)*sigmaP(:)-(sigmaXP)**2.d0)
    if(rmsEmittance(3) .ne. rmsEmittance(3)) then
        write(*,*) "Problem!", rmsEmittance(3), sigmaX(3), sigmaP(3), sigmaXP(3)
    end if
	write(41,*) time(i),j*nMacroParticle,rmsEmittance
end do
end subroutine emittance


!-------------------------------------------------------------------
! subroutine radialProfile
! Calculates radial profile of electron pulse
!-------------------------------------------------------------------
subroutine radialProfile(nStep,nProc)
implicit none
integer, intent(in) :: nStep, nProc

real, parameter :: sigmaXBin=10.0e-6 !m
integer :: nBin, nPart
integer :: i,j,stat, iBin
real :: x, xMin, xMax
real, allocatable :: radialData(:)
character(80) :: fileIn

xMin=1.0
xMax=0.0

do i=1,nProc
	stat=0
	if (i < 10) then
		write(fileIn,'(i,a,i1,a)') nStep,"-x-",i,".txt"
	else
		write(fileIn,'(i,a,i2,a)') nStep,"-x-",i,".txt"
	end if
	open(unit=40,file=trim(fileIn),status="old")
	do while(stat==0)
		read(40,*,iostat=stat) x
        if(stat == 0) then
            if(x > xMax) xMax=x
            if(x < xMin) xMin=x
        end if
	end do
	close(40)
end do

nBin=int((xMax-xMin)/sigmaXBin+1)
write(*,*) "N bins used",nBin
allocate(radialData(nBin))
radialData(:)=0
nPart=0
do j=1,nProc
    stat=0
    if (j < 10) then
        write(fileIn,'(i,a,i1,a)') nStep,"-x-",j,".txt"
    else
        write(fileIn,'(i,a,i2,a)') nStep,"-x-",j,".txt"
    end if
    open(unit=40,file=trim(fileIn),status="old")
    do while(stat==0)
        read(40,*,iostat=stat) x
        if(stat == 0) then
            nPart=nPart+1
            iBin=int((x-xMin)*nBin/(xMax-xMin))+1
            if(iBin > nBin) iBin=nBin
            radialData(iBin)=radialData(iBin)+1
        end if
    end do
    close(40)
end do

do i=1,nBin
    write(41,*) xMin+(xMax-xMin)/real(nBin)*i,radialData(i)/nPart
end do
end subroutine radialProfile

!-------------------------------------------------------------------
! subroutine longitudinalProfile
! Calculates radial profile of electron pulse
!-------------------------------------------------------------------
subroutine longitudinalProfile(nStep,nProc)
implicit none
integer, intent(in) :: nStep, nProc

real :: sigmaXBin=1.0e-6 !m
integer :: nBin, nPart
integer :: i,j,stat, iBin
real :: x, xMin, xMax, x1, x2, xAve
real, allocatable :: radialData(:)
character(80) :: fileIn

xMin=1.0
xMax=0.0
xAve = 0.0

do i=1,nProc
	stat=0
	if (i < 10) then
		write(fileIn,'(i,a,i1,a)') nStep,"-x-",i,".txt"
	else
		write(fileIn,'(i,a,i2,a)') nStep,"-x-",i,".txt"
	end if
	open(unit=40,file=trim(fileIn),status="old")
	do while(stat==0)
		read(40,*,iostat=stat) x1, x2, x
        if(stat == 0) then
            if(x > xMax) xMax=x
            if(x < xMin) xMin=x
        end if
	end do
	close(40)
end do

nBin=int((xMax-xMin)/sigmaXBin+1)
do while ( nBin > 1e4 )
    sigmaXBin = sigmaXBin * 10.0
    nBin=int((xMax-xMin)/sigmaXBin+1)
end do
do while ( nBin < 10 )
    sigmaXBin = sigmaXBin/10.0
    nBin=int((xMax-xMin)/sigmaXBin+1)
end do

write(*,*) "N bins used",nBin
allocate(radialData(nBin))
radialData(:)=0
nPart=0
do j=1,nProc
    stat=0
    if (j < 10) then
        write(fileIn,'(i,a,i1,a)') nStep,"-x-",j,".txt"
    else
        write(fileIn,'(i,a,i2,a)') nStep,"-x-",j,".txt"
    end if
    open(unit=40,file=trim(fileIn),status="old")
    do while(stat==0)
        read(40,*,iostat=stat) x1, x2, x
        if(stat == 0) then
            nPart=nPart+1
            iBin=int((x-xMin)*nBin/(xMax-xMin))+1
            if(iBin > nBin) iBin=nBin
            radialData(iBin)=radialData(iBin)+1
            xAve = xAve + x
        end if
    end do
    close(40)
end do
xAve = xAve/real(nPart)

do i=1,nBin
    write(41,*) xMin+(xMax-xMin)/real(nBin)*i-xAve,radialData(i)/nPart
end do
end subroutine longitudinalProfile

end program postprocessing
