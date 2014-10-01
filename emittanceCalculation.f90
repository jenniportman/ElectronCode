program emittance
implicit none
integer :: nSteps
character(10) :: stringField
integer :: i,j,stat,ierr
character(80) :: fileIn,fileOut,fileOut2
real, dimension(3) :: meanX,meanX2,meanP,meanP2,meanXP,rmsEmittance,xVect,pVect
real, dimension(3) :: sigmaX,sigmaP,sigmaXP
character(5), dimension(17) :: folderList
integer, dimension(17) :: mpList, neInit

write(*,*) "Insert total number of time steps"
read(*,*) nSteps
write(*,*) "Insert extraction field"
read(*,*) stringField

write(fileOut,*) "emittance_data_norm_",trim(stringField),".txt"
open(unit=41,file=trim(fileOut),status="unknown",recl=400)
write(fileOut2,*) "emittance_data_",trim(stringField),".txt"
open(unit=42,file=trim(fileOut2),status="unknown",recl=400)

folderList=(/ "1e5","2e5","5e5", "1e6", "2e6", "3e6", "5e6", "7e6", "1e7", "1.5e7", "2e7", "3e7", "4e7", "6e7", "8e7", "1e8", "2e8" /)
mpList=(/ (10,i=1,3),(100,i=1,14) /)
neInit=(/ 1e5, 2e5, 5e5, 1e6, 2e6, 3e6, 5e6, 7e6, 1e7, 1.5e7, 2e7, 3e7, 4e7, 6e7, 8e7, 1e8, 2e8 /)

do i=1,size(folderList)
    meanX=0.0
    meanX2=0.0
    meanP=0.0
    meanP2=0.0
    meanXP=0.0
	j=0
	stat=0
	write(fileIn,'(a,a,i3,a)') trim(folderList(i)),"/",nSteps,"-x-all.txt"
    ierr=0
	open(unit=40,file=trim(fileIn),status="old",iostat=ierr)
    if(ierr == 0) then
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
        write(*,*) "Nel", j
        meanP(:)=meanP(:)/mpList(i)
        meanP2(:)=meanP2(:)/mpList(i)**2.d0
        meanXP(:)=meanXP(:)/mpList(i)
        sigmaX(:)=(real(j)*meanX2(:)-meanX(:)*meanX(:))/(real(j)*real(j))
        sigmaP(:)=(real(j)*meanP2(:)-meanP(:)*meanP(:))/(real(j)*real(j))
        sigmaXP(:)=(real(j)*meanXP(:)-meanX(:)*meanP(:))/(real(j)*real(j))
        rmsEmittance(:)=sqrt(sigmaX(:)*sigmaP(:)-(sigmaXP)**2.d0)
        if(rmsEmittance(3) .ne. rmsEmittance(3)) then
            write(*,*) "Problem!", rmsEmittance(3)
        end if
        write(41,*) neInit(i),j*mpList(i),rmsEmittance, sqrt(sigmaX(:))
        write(42,*) neInit(i),j*mpList(i),rmsEmittance(:)/meanP(3)*real(j), sqrt(sigmaX(:))
    end if
end do

close(41)
close(42)

end program emittance
