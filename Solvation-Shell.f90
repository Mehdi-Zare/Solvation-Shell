program fun
implicit none

integer  ,  parameter                           ::      dp = selected_real_kind(15, 307)
character(len=50)                               ::      file1 
character(len=50)                               ::      file2 
character(len=50)                               ::      filename
integer                                         ::      error_flag, alloc_err                   
integer                                         ::      iii,jjj,kkk,nline1, nline2, totlines
integer,   parameter                            ::      metalMM=1229
integer,   parameter                            ::      metalQM=51
integer,   parameter                            ::      Ads=10
integer,   parameter                            ::      Water=2200*3
real(dp), dimension(2)                          ::      Xox, Yox, Zox  ! Oxygen atom coordinates            
real(dp), dimension(6600)                       ::      Xwi,Ywi,Zwi      ! i is for initial and f for final
real(dp), dimension(10)                         ::      Xai,Yai,Zai
character(len=256), dimension(6600)             ::      waterID
character(len=256), dimension(10)               ::      AdsID
real(dp)                                        ::      a               ! for FLAG
real                                            ::      HZ, LZ, IZ
integer                                         ::      nBin,totalimage
real(dp),allocatable, dimension(:)              ::      Bins, SumBins
real(dp)                                        ::      lowerbond, upperbond, Volume

totalimage=1000
HZ=15
LZ=0
IZ=0.1
nBin=((HZ-LZ)/IZ) ! the first 1 because fortran round the number down, the
                        !second 1 for the bin above HZ
!write(*,*) "nBin=" , nBin
allocate(Bins(nBin))
allocate(SumBins(nBin))
SumBins=0.0_dp

do iii=1,totalimage
  write(*,*) iii

  write(file1,100) iii
  100 format('image-', I5.5)
    call countline(file1,nline1)          
         if ( error_flag == 1  ) then
            call EXIT(0)           !FLAG
         end if
    totlines=metalMM+Water+metalQM+Ads+2
         if ( totlines /= nline1 ) then
                  write(*,*) "Total lines in ", file1, " must be  ", totlines
                  call EXIT(0)    !FLAG
         end if 
  call ReadConfig(file1,metalMM,metalQM,Ads,Water,waterID,Xwi,Ywi,Zwi,AdsID,Xai,Yai,Zai)
  call RotationalCorrelationFunction(Ads,Water,waterID,Xwi,Ywi,Zwi,AdsID,Xai,Yai,Zai,HZ,LZ,IZ,nBin,Bins)
  

  lowerbond=LZ
  upperbond=lowerbond+IZ
 ! Volume=(((4.0_dp/3.0_dp)*(3.14_dp))*((upperbond**3)-(lowerbond**3)))/2.0_dp ! b/c half is empty of watre(just metal atoms) 
  do jjj=1,nBin
    !SumBins(jjj)=SumBins(jjj)+(Bins(jjj)/Volume)
    SumBins(jjj)=SumBins(jjj)+Bins(jjj)
  end do  

end do

open    ( unit = 5000, file = 'RESULTS', status = 'new')
write   ( 5000, 5010)
5010 format ( 3x, "Bin No.", 11x, "Density(#H2O/Ang^3)", 16x,"LowerBond", 16x, "UpperBond", 16x, "Volume")
write   (5000,5020)
5020 format (2x, "============================================================",&
                "==============================================================")
lowerbond=LZ
upperbond=lowerbond+IZ
Volume=(((4.0_dp/3.0_dp)*(3.14_dp))*((upperbond**3)-(lowerbond**3)))/2.0_dp !b/c half is empty of watre(just metal atoms)
do jjj=1,nBin
    write(5000,5030) jjj , SumBins(jjj)/REAL(totalimage)/Volume, lowerbond, upperbond, Volume
    5030 format (3x, I3, 15x ,ES20.13 , 16x, F6.2, 17x, F6.2, 17x, ES20.13)
    lowerbond=upperbond
    upperbond=lowerbond+IZ
    Volume=(((4.0_dp/3.0_dp)*(3.14_dp))*((upperbond**3)-(lowerbond**3)))/2.0_dp ! b/c half is empty of watre(just metal atoms)
end do
close(5000)

deallocate(Bins,SumBins, stat =alloc_err)


end program fun


subroutine RotationalCorrelationFunction(Ads,Water,waterID,Xwi,Ywi,Zwi,AdsID,Xai,Yai,Zai,HZ,LZ,IZ,nBin,Bins)
implicit none

integer  ,  parameter                           ::      dp = selected_real_kind(15, 307)

!local variable
integer                                         ::      ierror,  error_flag, alloc_err
integer                                         ::      ii, jj, kk
real(dp), dimension(2)                          ::      Xox, Yox, Zox  ! Oxygenatom coordinates
real(dp)                                        ::      distance1i, distance2i, distance1f, distance2f, unitdipol
real(dp)                                        ::      dotdipolXi, dotdipolYi, dotdipolZi, dotdipolXf
real(dp)                                        ::      dotdipolYf,dotdipolZf,lengthi,lengthf
real(dp)                                        ::      Xi,Yi,Zi, Xf, Yf, Zf
real(dp), dimension(Water)                      ::      dipoli, dipolf
real                                            ::      lowerbond, upperbond
!integer,   parameter                            ::      metalMM=1229
!integer,   parameter                            ::      metalQM=51
!integer,   parameter                            ::      Ads=10
!integer,   parameter                            ::      Water=2200

!parameter types and definition
integer,                              intent(in)                  :: Ads,Water, nBin
real(dp), dimension(Water),           intent(in)                  :: Xwi,Ywi,Zwi
character(len=256),dimension(Water),  intent(in)                  :: waterID
real(dp), dimension(Ads),             intent(in)                  :: Xai,Yai,Zai
character(len=256),dimension(Ads),    intent(in)                  :: AdsID
real,                                 intent(in)                  :: HZ, LZ, IZ
real(dp) , dimension(nBin),            intent(out)                 :: Bins

! Get the Oxygen atom coordinate of the adsorbate, it is the same for all
! conformations

jj=1 
do ii=1,Ads
        if (AdsID(ii) == 'O' ) then
        Xox(jj)=Xai(ii); Yox(jj)=Yai(ii); Zox(jj)=Zai(ii)
        jj=jj+1
        end if
end do


!allocate(Bins(nBin))
Bins=0.0_dp

!write(*,*) Xox, Yox, Zox  !FLAG

!dipoli=0.0_dp
!dipolf=0.0_dp
kk=1 ! controls the number of water molecule in 5 ang at time zero (initial)
jj=1 ! control the number of water molecules that still exist in 5 ang at time t (final)

lowerbond=LZ
upperbond=lowerbond+IZ
do jj=1,nBin
  do ii=1,Water
        if (waterID(ii) == 'Ow') then
                distance1i=SQRT(((Xwi(ii)-Xox(1))**2) + ((Ywi(ii)-Yox(1))**2)  + ((Zwi(ii)-Zox(1))**2) )
                distance2i=SQRT(((Xwi(ii)-Xox(2))**2) + ((Ywi(ii)-Yox(2))**2)  + ((Zwi(ii)-Zox(2))**2) )
                if ( (distance1i >= lowerbond .and. distance1i < upperbond)&
                 .or. (distance2i >= lowerbond .and. distance2i < upperbond) ) then
                        Bins(jj)=Bins(jj)+1
                end if
        end if        
  end do
  lowerbond=upperbond
  upperbond=lowerbond+IZ

end do

WRITE(*,*)'Length of Bins vecotr is',  SHAPE(Bins)
end subroutine RotationalCorrelationFunction

!deallocate(waterID,Xw,Yw,Zw,AdsID,Xa,Ya,Za,stat = alloc_err)

!!!!!!!!!!!!!!!!!!!!!!!!!!!! SUBROUTINE FOR READING CONFIG FORMAT FILE !!!!!!!!!!!!!!!!!!!!!!
subroutine ReadConfig(filename,metalMM,metalQM,Ads,Water,waterID,Xwat,Ywat,Zwat,AdsID,Xads,Yads,Zads)
implicit none

integer  ,  parameter                           ::      dp = selected_real_kind(15, 307)

!local variable
integer                                         ::      ierror,  error_flag, alloc_err
integer                                         ::      ii
character(len=256)                              ::      skip ! The variable forskipping lines
character(len=256)                              ::      str_2, str_3, str_4 !The variables for skipping strings in a line
!integer,   parameter                            ::      metalMM=1229
!integer,   parameter                            ::      metalQM=51
!integer,   parameter                            ::      Ads=10
!integer,   parameter                            ::      Water=2200

!parameter types and definition
character(len=30)    ,                intent(in)                  ::      filename 
integer,                              intent(in)                  ::      metalMM,metalQM,Ads,Water
real(dp), dimension(Water),  intent(out)                          ::      Xwat,Ywat,Zwat
character(len=256),dimension(Water), intent(out)                  ::      waterID
real(dp), dimension(Ads),  intent(out)                            ::      Xads,Yads,Zads 
character(len=256),dimension(Ads), intent(out)                    ::      AdsID


!Read the energies in the file
!allocate(metalMMID(metalMM))
!allocate(XmMM(metalMM))
!allocate(YmMM(metalMM))
!allocate(ZmMM(metalMM))
!allocate(waterID(Water))
!allocate(Xwat(waterlines))
!allocate(Ywat(waterlines))
!allocate(Zwat(waterlines))

!allocate(metalQMID(metalQM))
!allocate(XmQM(metalQM))
!allocate(YmQM(metalQM))
!allocate(ZmQM(metalQM))

!allocate(AdsID(Ads))
!allocate(Xads(Ads))
!allocate(Yads(Ads))
!allocate(Zads(Ads))

error_flag = 0

open(unit=77777, file=filename, status='old', action='read', iostat = ierror)

! skip the header lines, 2 lines
do ii=1,2
   read(77777,'(A)', iostat = ierror) skip
        if (ierror /=0 ) then                  !!!!FALG
          write(*,*) " Problem with reading file: ", filename
          error_flag = 1
          call EXIT(0)
        end if
  ! write(*,*) skip
end do
! Read Metal MM atoms, those are 1229 Pt atoms
do ii=1,metalMM
        read(77777,*, iostat = ierror) skip
                if (ierror /=0 ) then                  !!!!FALG
                write(*,*) " Problem with reading file: ", filename
                error_flag = 1
                call EXIT(0)
                end if
        !write(*,*)  skip
end do

do ii=1,Water
        read(77777,*, iostat = ierror) waterID(ii), str_2, str_3, str_4
                if (ierror /=0 ) then                  !!!!FALG
                write(*,*) " Problem with reading file: ", filename
                error_flag = 1
                call EXIT(0)
                end if
!        write(*,*)  waterID(ii)
        read(str_2,*) Xwat(ii)
!        write(*,*) Xwat(ii)
        read(str_3,*) Ywat(ii)
!        write(*,*) Ywat(ii)
        read(str_4,*) Zwat(ii)
!        write(*,*) Zwat(ii)
end do

do ii=1,metalQM
        read(77777,*, iostat = ierror) skip
                if (ierror /=0 ) then                  !!!!FALG
                write(*,*) " Problem with reading file: ", filename
                error_flag = 1
                call EXIT(0)
                end if
        !write(*,*) skip
end do
do ii=1,Ads
        read(77777,*, iostat = ierror) AdsID(ii), str_2, str_3, str_4
                if (ierror /=0 ) then                  !!!!FALG
                write(*,*) " Problem with reading file: ", filename
                error_flag = 1
                call EXIT(0)
                end if
!        write(*,*)  AdsID(ii)
        read(str_2,*) Xads(ii)
!        write(*,*) Xads(ii)
        read(str_3,*) Yads(ii)
!        write(*,*) Yads(ii)
        read(str_4,*) Zads(ii)
!        write(*,*) Zads(ii)
end do


close(77777)

!deallocate(metalMMID, XmMM, YmMM, ZmMM,waterID, stat = alloc_err)
!deallocate(metalQMID, XmQM, YmQM, ZmQM,AdsID, stat = alloc_err)

end subroutine ReadConfig
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! aEND OF SUBROUTINE energyValues !!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SUBROUTNE FOR COUNTING lines IN A FILE!!!!!!!!!!!!!
subroutine countline(filename, nline)
implicit none

integer  ,  parameter                   ::      dp = selected_real_kind(15, 307)

!local variable
integer                                 ::      ierror,  error_flag

!parameter types and definition
character(len=30)    ,  intent(in)                  ::      filename
integer              ,  intent(out)                 ::      nline


error_flag = 0

!read the # of lines in the file 
nline=0
open(1000, file=filename, status='old', action='read', iostat = ierror)
do
 read (1000, *, end=10)
        if (ierror /=0 ) then                  !!!!FALG
          write(*,*) " Problem with counting lines in file:  " , filename
           error_flag = 1
          call EXIT(0)
        end if
 nline=nline+1
end do
10 close(1000)

end subroutine countline
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! END OF SUBROUTINE COUNTLINE !!!!!!!!!!!!!!!!!!!!!

