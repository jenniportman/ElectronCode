program generate
implicit none
real, parameter :: pi=3.1415926535897932384626,c=2.99792458E8
integer :: nElectrons, nParts, nMPR, nProfile
real :: tSigma, xSigma, Ef, W, Eph
character(80) :: fileIn, fileSlice, fileOut 
character(10) :: stringNElectrons,stringNMPR, fileNameAppend

integer :: i,j,nEGen,nCheck, ierr, iParts
real, allocatable :: probability(:)
real :: dt, x, y, z, px, py, pz, Eel, theta, phi, thetaMax, massEl, vel

integer :: iClock
integer,dimension(8) :: iSeed
real :: rand,urand,vrand,wrand
real :: r, angle

call system_clock(iClock)
iSeed(:)=iClock
call random_seed(PUT=iseed)
call random_number(rand)

write(fileIn,*) "input.txt"
write(fileSlice,*) "probability.txt"

write(fileNameAppend,*) " "
open(40,file=fileIn,status="old")
read(40,*) nElectrons
read(40,*) nParts
read(40,*) nMPR
read(40,*) tSigma
read(40,*) xSigma
read(40,*) Ef,W,Eph
read(40,*) nProfile
read(40,*,iostat=ierr) fileNameAppend
close(40)

dt=3.0*tSigma/(nParts/2)
massEl=0.510998928

! read in file with probabilities for number of electron for each time slice
allocate(probability(nParts))
open(41,file=fileSlice,status="old")
do i=1,nParts
	read(41,*) probability(i)
end do
close(41)

write(stringNElectrons,'(i10)') nElectrons
write(stringNMPR,'(i5)') nMPR
if(nProfile == 0) then
    if( fileNameAppend == "") then
        write(fileOut,*) "IniCondition-G-",trim(adjustl(stringNElectrons)),"-",trim(adjustl(stringNMPR)),".txt"
    else
        write(fileOut,*) "IniCondition-G-",trim(adjustl(stringNElectrons)),"-",trim(adjustl(stringNMPR)),"-",trim(adjustl(fileNameAppend)),".txt"
    end if
else if(nProfile == 1) then
    if( fileNameAppend == '') then
        write(fileOut,*) "IniCondition-U-",trim(adjustl(stringNElectrons)),"-",trim(adjustl(stringNMPR)),".txt"
    else
        write(fileOut,*) "IniCondition-U-",trim(adjustl(stringNElectrons)),"-",trim(adjustl(stringNMPR)),"-",trim(adjustl(fileNameAppend)),".txt"
    end if
else if(nProfile == 2) then
    if( fileNameAppend == '') then
        write(fileOut,*) "IniCondition-E-",trim(adjustl(stringNElectrons)),"-",trim(adjustl(stringNMPR)),".txt"
    else
        write(fileOut,*) "IniCondition-E-",trim(adjustl(stringNElectrons)),"-",trim(adjustl(stringNMPR)),"-",trim(adjustl(fileNameAppend)),".txt"
    end if
else
    write(*,*) "Error in nProfile"
    stop
end if

open(42,file=fileOut,status="unknown",recl=400)
nCheck=0
iParts=1
do i=1,nParts
	nEGen=nint(probability(i)*nElectrons)
	nCheck=nCheck+nEGen
    if(nEGen > 0) iParts=iParts+1
	do j=1,nEGen

		! generate momenta using 3 step model of photoemission
		call random_number(rand)
		Eel=rand*(Eph-W)+Ef+W
		thetaMax=acos(sqrt((Ef+W)/Eel))
		call random_number(theta)
		theta=theta*thetaMax
		call random_number(phi)
		phi=phi*2.0*pi
		vel=sqrt(2*Eel/massEl)
		px=massEl*nMPR*vel*sin(theta)*cos(phi)
		py=massEl*nMPR*vel*sin(theta)*sin(phi)
		pz=vel*cos(theta)
        if(pz < sqrt(2*(Ef+W)/massEl)) write(*,*) "ERROR",pz,sqrt(2*(Ef+W)/massEl)
		pz=massEl*nMPR*sqrt(pz**2.0-2.0*(Ef+W)/massEl)

		if(nProfile == 0) then
			! generate x,y positions with normal distribution using Box-Muller transform
			wrand=0.0
			do while(wrand==0 .or. wrand>=1)
				call random_number(urand)
				urand=urand*2.0-1.0
				call random_number(vrand)
				vrand=vrand*2.0-1.0
				wrand=urand**2.0+vrand**2.0
			end do
			wrand=sqrt(-2.0*log(wrand)/wrand)
			x=urand*wrand*xSigma
			y=vrand*wrand*xSigma
            z=0
        else if(nProfile == 1) then
			! generate x,y positions with uniform distribution
			call random_number(urand)
			call random_number(vrand)
            angle = urand*2*pi
            ! numerical factors ensure that we get standard deviation xSigma
            ! sqrt vrand comes from generating points uniformly in circle
            r = sqrt(vrand)*xSigma*2.0
			x = r * cos(angle)
			y = r * sin(angle)
            z=0
        else if(nProfile == 2) then
			call random_number(urand)
			call random_number(vrand)
            ! numerical factors ensure that we get standard deviation xSigma
            ! expression inside sqrt comes from generating r with correct distribution
            r = xSigma * sqrt(5.0) * sqrt(1.0-(1.0-urand)**(2.0/3.0))
            angle = vrand*2*pi
			x = r * cos(angle)
			y = r * sin(angle)
            z=0
		end if

		write(42,*) (i-1)*dt,dt,nEGen,x,y,z,px,py,pz
	end do
end do
write(*,*) nCheck
write(*,*) "electron generation steps",iParts-1
write(42,*) "The above data file has the following format:"
write(42,*) "time, time step size, number of particle generated in this time step, x, y, z, px, py, pz"
write(42,*) "total particles generated:",nCheck
write(42,*) "total time steps",nParts, "actual time steps",iParts-1
write(42,*) "the following parameters were used"
write(42,*) "dx(m), dy(m)",xSigma, xSigma
write(42,*) "Ef(MeV),Eph(MeV),W(MeV)",Ef,Eph,W
write(42,*) "mpr=",nMPR
write(42,*) "me(MeV/c^2)",massEl
close(42)

end program generate
