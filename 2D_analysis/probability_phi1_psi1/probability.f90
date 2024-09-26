PROGRAM mtd_probability
IMPLICIT NONE
REAL*8 :: gridmin1, gridmax1, griddiff1, gridmin2, gridmax2, griddiff2
REAL*8 :: gridmin3, gridmax3, griddiff3, gridmin4, gridmax4, griddiff4
REAL*8 :: kt0, kt, ktb, bias_fact, alpha_cv, deltaT, den
REAL*8 :: diff_s2, ds2, ss, hh, dum, num, gamma_,s1,s2
REAL*8, ALLOCATABLE :: prob(:,:)
REAL*8, ALLOCATABLE :: cv1(:), cv2(:), vbias(:)
REAL*8, ALLOCATABLE :: cv3(:), cv4(:), ct(:)
INTEGER :: mtd_steps, md_steps, i, t_min, t_max
INTEGER :: mtd_max, w_cv, w_hill, i_mtd, i_md, s3, s4
INTEGER :: i_s1, i_s2, i_s3, i_s4, nbin1, nbin2, nbin3, nbin4
INTEGER :: index1, index2, index3, index4

REAL*8, PARAMETER :: kb=1.9872041E-3 !kcal K-1 mol-1
REAL*8, PARAMETER :: kj_to_kcal = 0.239006

OPEN(1,FILE='input',STATUS='old')
!OPEN(2,FILE='PROB.dat',STATUS='replace',form='unformatted')
OPEN(3,FILE='Pu.dat',STATUS='unknown')
OPEN(11,FILE='COLVAR',STATUS='old')
OPEN(12,FILE='HILLS',STATUS='old')
OPEN(21,file='data_ct.dat',status='old')
OPEN(22,file='data_vbias.dat',status='old')


CALL get_steps(11,md_steps)
CALL get_steps(12,mtd_steps)

print *, 'md_steps=', md_steps, 'mtd_steps=', mtd_steps

read(1,*) kt0, kt, bias_fact
read(1,*) t_min, t_max

IF(t_max.gt.md_steps)STOP '!!ERROR: t_max > total MD steps'

read(1,*) gridmin1, gridmax1, griddiff1
read(1,*) gridmin2, gridmax2, griddiff2
read(1,*) w_cv, w_hill

! deltaT=(bias_fact - 1.d0)*kt
! deltaT = (bias_fact - 1.d0)*kt0
! alpha_cv = (kt + deltaT)/deltaT

kt = kb*kt
! gamma_ = 1/alpha
gamma_ = (bias_fact - 1.0)/bias_fact
write(*,*) 'gamma_=', gamma_


WRITE(*,'(A,I10)')'No: of MTD steps        =',mtd_steps
WRITE(*,'(A,I10)')'No: of MD  steps        =',md_steps
WRITE(*,'(A,F9.2)')'Physical Temp (K)      =',kt0
WRITE(*,'(A,F9.2)')'CV Temp (K)            =',kt
WRITE(*,'(A,F9.2)')'Bias Factor (K)        =',bias_fact
WRITE(*,'(A,I10)')'Print Freq. cvmdck_mtd  =',w_cv
WRITE(*,'(A,I10)')'Freq. of Hill Update    =',w_hill
WRITE(*,'(A,I10)')'Reweigtht: Min step     =',t_min
WRITE(*,'(A,I10)')'Reweigtht: Max step     =',t_max


ALLOCATE(cv1(md_steps),cv2(md_steps))
ALLOCATE(vbias(md_steps),ct(mtd_steps))


DO i_md=1,md_steps
 READ(11,*) dum, cv1(i_md),cv2(i_md),dum,dum,dum,dum,dum,dum 
     IF( cv1(i_md) .gt.  3.14d0)  cv1(i_md) = cv1(i_md) - 6.28d0
     IF( cv1(i_md) .lt. -3.14d0 ) cv1(i_md) = cv1(i_md) + 6.28d0
     IF( cv2(i_md) .gt.  3.14d0)  cv2(i_md) = cv2(i_md) - 6.28d0
     IF( cv2(i_md) .lt. -3.14d0 ) cv2(i_md) = cv2(i_md) + 6.28d0
END DO

nbin1 = NINT((gridmax1-gridmin1)/griddiff1)+1
nbin2 = NINT((gridmax2-gridmin2)/griddiff2)+1

write(*,*) nbin1, nbin2

write(*,*) "reading vbias.dat file"

do i_md=1,md_steps
  read(22,*) dum, vbias(i_md)
end do

write(*,*) "reading ct.dat file"
do i_mtd=1,mtd_steps
  read(21,*) dum, ct(i_mtd)
end do

ALLOCATE(prob(nbin1,nbin2))


WRITE(*,*) 'calculating  probability'

den=0.d0
prob=0.d0

DO i_md=1,md_steps
 IF((i_md.GT.t_min).AND.(i_md.LT.t_max))THEN
 index1 = nint((cv1(i_md)-gridmin1)/griddiff1) +1
 index2 = nint((cv2(i_md)-gridmin2)/griddiff2) +1

     IF(index1.gt.0.and.index2.gt.0.and.index1.le.nbin1.and.index2.le.nbin2) then

    ! i_mtd=(i_md*w_cv/w_hill)+1
      i_mtd=((i_md-1)*w_cv/w_hill)+1
      dum=vbias(i_md) - ct(i_mtd)
      prob(index1,index2) = prob(index1,index2) + DEXP(dum/kt)
      den=den+DEXP(dum/kt)
     END IF

 END IF
END DO

  den=den*griddiff1*griddiff2


DO i_s1=1,nbin1
   s1=DFLOAT(i_s1-1)*griddiff1+gridmin1
   DO i_s2=1,nbin2
      s2=DFLOAT(i_s2-1)*griddiff2+gridmin2
      prob(i_s1,i_s2)=prob(i_s1,i_s2)/den
          
    ! WRITE(2)prob(i_s1,i_s2)
      WRITE(3,'(3E16.8)')s1,s2,prob(i_s1,i_s2)
   END DO
   WRITE(3,*)
END DO
      WRITE(*,'(A)')'Unbiased distribution written in Pu.dat'

close(1)
close(2)
close(3)
close(11)
close(12)
close(21)
close(22)

DEALLOCATE(cv1, cv2, vbias, ct)
DEALLOCATE(prob)



END PROGRAM

SUBROUTINE get_steps(iunit,nsteps)
IMPLICIT NONE
INTEGER :: iunit, nsteps
INTEGER :: ios
nsteps=0
REWIND(iunit)
  Read_Loop: DO
  READ(iunit,*,IOSTAT=ios)
    IF(ios.ne.0)EXIT Read_Loop
    nsteps=nsteps+1
  END DO Read_Loop
REWIND(iunit)
END SUBROUTINE


