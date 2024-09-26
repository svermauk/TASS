! MPi version of fortran code to calculate vbaias from TASS data
!
! mpif90 vbias_factor_mpi.f90 -o vbias_mpi.x
!
! mpirun -np 2 ./vbias_mpi.x
!
!
! Modified by Anji Babu Kapakayala and Shivani Verma
!
PROGRAM vbias_factor
IMPLICIT NONE
REAL*8 :: gridmin1, gridmax1, griddiff1, gridmin2, gridmax2, griddiff2
REAL*8 :: gridmin3, gridmax3, griddiff3, gridmin4, gridmax4, griddiff4
REAL*8 :: kt0, kt, ktb, bias_fact, alpha_cv, deltaT, den
REAL*8 :: diff_s2, ds2, ss, hh, dum, dum_vb, num, gamma_
REAL*8, ALLOCATABLE :: cv1(:), cv2(:), ht(:), vbias(:), hill(:)
REAL*8, ALLOCATABLE :: cv3(:), cv4(:), width(:)
INTEGER :: mtd_steps, md_steps, i, t_min, t_max 
INTEGER :: mtd_min, mtd_max, w_cv, w_hill, i_mtd, i_md, nbin1
INTEGER :: i_s1, i_s2, i_s3, i_s4, nbin2, nbin3, nbin4

REAL*8, PARAMETER :: kb=1.9872041E-3 !kcal K-1 mol-1
REAL*8, PARAMETER :: kj_to_kcal = 0.239006

! MPI variables
logical :: parent
integer :: rank, icpu, ncpu,a

CALL MPI_Start
CALL Set_Parent(parent)

if (parent) then

open(1,FILE='input',STATUS='old')
OPEN(11,FILE='COLVAR',STATUS='old')
OPEN(12,FILE='HILLS',STATUS='old')
OPEN(21,FILE='data_vbias.dat',STATUS='replace')

CALL get_steps(11,md_steps)
!

CALL get_steps(12,mtd_steps)
!
print *, 'md_steps=', md_steps, 'mtd_steps=', mtd_steps

read(1,*) kt0, kt, bias_fact
read(1,*) t_min, t_max

IF(t_max.gt.md_steps)STOP '!!ERROR: t_max > total MD steps'
md_steps=t_max

read(1,*) gridmin1, gridmax1, griddiff1
read(1,*) gridmin2, gridmax2, griddiff2
read(1,*) w_cv, w_hill

! deltaT=(bias_fact - 1.d0)*kt
! deltaT = (bias_fact - 1.d0)*kt0
! alpha_cv = (kt + deltaT)/deltaT
endif

!MPI routines
!broadcast ncv
!CALL IBcast_1(ns,1)
CALL IBcast_1(w_cv,1)
CALL IBcast_1(w_hill,1)
CALL IBcast_1(t_min,1)
CALL IBcast_1(t_max,1)
CALL IBcast_1(mtd_steps,1)
CALL IBcast_1(md_steps,1)
CALL RBcast_1(kt0,1)
CALL RBcast_1(kt,1)
CALL RBcast_1(bias_fact,1)



if (parent) then

WRITE(*,'(A,I10)')'No: of MTD steps        =',mtd_steps
WRITE(*,'(A,I10)')'No: of MD  steps        =',md_steps
WRITE(*,'(A,F9.2)')'Physical Temp (K)      =',kt0
WRITE(*,'(A,F9.2)')'CV Temp (K)            =',kt
WRITE(*,'(A,F9.2)')'Bias Factor (K)        =',bias_fact
WRITE(*,'(A,I10)')'Print Freq. cvmdck_mtd  =',w_cv
WRITE(*,'(A,I10)')'Freq. of Hill Update    =',w_hill
WRITE(*,'(A,I10)')'Reweigtht: Min step     =',t_min
WRITE(*,'(A,I10)')'Reweigtht: Max step     =',t_max

endif

kt = kb*kt
! gamma_ = 1/alpha
gamma_ = (bias_fact - 1.0)/bias_fact
if (parent) write(*,*) 'gamma_=', gamma_


ALLOCATE(cv1(md_steps),cv2(md_steps))
ALLOCATE(vbias(md_steps))
ALLOCATE(ht(mtd_steps))
ALLOCATE(hill(mtd_steps),width(mtd_steps))


if (parent) then
DO i_md=1,md_steps
 READ(11,*) dum, cv1(i_md),dum,cv2(i_md),dum,dum,dum,dum,dum 
     IF( cv1(i_md) .gt.  3.14d0)  cv1(i_md) = cv1(i_md) - 6.28d0
     IF( cv1(i_md) .lt. -3.14d0 ) cv1(i_md) = cv1(i_md) + 6.28d0
     IF( cv2(i_md) .gt.  3.14d0)  cv2(i_md) = cv2(i_md) - 6.28d0
     IF( cv2(i_md) .lt. -3.14d0 ) cv2(i_md) = cv2(i_md) + 6.28d0
END DO


nbin1 = NINT((gridmax1-gridmin1)/griddiff1)+1
nbin2 = NINT((gridmax2-gridmin2)/griddiff2)+1

write(*,*) nbin1, nbin2

DO i_mtd=1,mtd_steps
 READ(12,*) dum,hill(i_mtd),width(i_mtd),ht(i_mtd)
      IF( hill(i_mtd) .gt.  3.14d0) hill(i_mtd) = hill(i_mtd) - 6.28d0
      IF( hill(i_mtd) .lt. -3.14d0 )hill(i_mtd) = hill(i_mtd) + 6.28d0
       ht(i_mtd)=ht(i_mtd)*kj_to_kcal
END DO


 write(*,*) 'calculating vbias'
endif

! Broadcasting Grid data and bin info
CALL RBcast_1(gridmin1,1)
CALL RBcast_1(gridmax1,1)
CALL RBcast_1(griddiff1,1)
CALL RBcast_1(gridmin2,1)
CALL RBcast_1(gridmax2,1)
CALL RBcast_1(griddiff2,1)
CALL IBcast_1(nbin1,1)
CALL IBcast_1(nbin2,1)

! Broadcasting HILL data
CALL RBcast(hill,mtd_steps)
CALL RBcast(width,mtd_steps)
CALL RBcast(ht,mtd_steps)

CALL Get_ncpu(ncpu)
CALL Get_cpuid(icpu)

if (parent)  write(*,*)'NCPUS',ncpu

! Broadcasting CV data
CALL RBcast(cv1,md_steps)
CALL RBcast(cv2,md_steps)
CALL RBcast_1(gamma_,1)
CALL RBcast_1(dum_vb,1)
!CALL RBcast(grid,ns*nbin(1))

!DO i_md=t_min,t_max
DO i_md=1,md_steps
!i_md=10202

  CALL  Distribute_Mtd_Min_Max(i_md,md_steps,w_cv,w_hill,mtd_min,mtd_max)

! print*,"cpu id:",icpu,"i_md",i_md,"mtd_min:",mtd_min,"mtd_max:",mtd_max

! i_md = 1_md -1
!td_max=((i_md-1)*w_cv/w_hill)+1
!F( MOD((i_md-1)*w_cv,w_hill).eq.0)  mtd_max=mtd_max-1
 ss=cv2(i_md)
 dum=0.d0
!-----------!  
 DO i_mtd=mtd_min,mtd_max
    ds2=width(i_mtd)*width(i_mtd)
    hh=ht(i_mtd)*gamma_
    diff_s2=ss-hill(i_mtd)
     if (diff_s2 .gt. 3.14d0 ) diff_s2 =diff_s2 - 6.28d0
     if (diff_s2 .lt.-3.14d0 ) diff_s2 =diff_s2 + 6.28d0
    diff_s2=diff_s2*diff_s2*0.5D0
    diff_s2=diff_s2/ds2
!   IF(ABS(diff_s2).GT.1.e-12)dum=dum+hh*DEXP(-diff_s2)
     dum=dum+hh*DEXP(-diff_s2)
!     print*, "i_mtd =",i_mtd ,"dum =", dum
!print*,"i_md:",i_md,"i_mtd:",i_mtd, "CPU:",icpu,"dum=:",dum
   END DO
!print*, "CPU:",icpu,"dum=:",dum
!------------!
  call Sync_Procs
  call GlobSumR_1(dum,dum_vb,1)

  if(parent) then 
!    vbias(i_md)=dum_vb
!    write(21,*) i_md, vbias(i_md)
    write(21,*) i_md, dum_vb
!   i_md = i_md + 1 
  endif

END DO
!100 FORMAT(I10,2X,F10.4,2X)
call Sync_Procs

if(parent) then
write(*,*) " vbias written in data_vbias.dat file. "
endif
close(1)
close(11)
close(21)
close(12)

DEALLOCATE(cv2, ht, vbias, hill)
DEALLOCATE( width, cv1)
call MPI_Stop

END PROGRAM
!--------------------------------------------------------------------!

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

!****************************************************************************************!
!   MPI ROUTINES
!****************************************************************************************!

SUBROUTINE MPI_Start()
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: i_err
!#if defined (_PARALLEL)
  call MPI_INIT(i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE MPI_Stop()
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: i_err
!#if defined (_PARALLEL)
call MPI_FINALIZE(i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE get_ncpu(ncpu)
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: ncpu, i_err
!ncpu=1
!#if defined (_PARALLEL)
call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE get_cpuid(icpu)
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: icpu, i_err
!icpu=0
!#if defined (_PARALLEL)
call MPI_COMM_RANK(MPI_COMM_WORLD,icpu,i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE IBcast(myint,leng)
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: leng, myint(*), i_err
!#if defined (_PARALLEL)
CALL MPI_BCAST(myint,leng,MPI_INTEGER,0,MPI_COMM_WORLD,i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE IBcast_1(myint,leng)
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: leng, myint, i_err
!#if defined (_PARALLEL)
CALL MPI_BCAST(myint,leng,MPI_INTEGER,0,MPI_COMM_WORLD,i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE RBcast(myreal,leng)
IMPLICIT NONE
INCLUDE 'mpif.h'
!INTEGER :: myreal(*)
REAL*8 :: myreal(*)
INTEGER :: leng, i_err
!#if defined (_PARALLEL)
CALL MPI_BCAST(myreal,leng,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE RBcast_1(myreal,leng)
IMPLICIT NONE
INCLUDE 'mpif.h'
!INTEGER :: myreal(*)
REAL*8 :: myreal
INTEGER :: leng, i_err
!#if defined (_PARALLEL)
CALL MPI_BCAST(myreal,leng,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE Sync_Procs
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER i_err
!#if defined (_PARALLEL)
call MPI_Barrier(MPI_COMM_WORLD,i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE Set_Parent(parent)
IMPLICIT NONE
INCLUDE 'mpif.h'
LOGICAL :: parent
INTEGER :: icpu, i_err
parent=.false.
!icpu=0
!#if defined (_PARALLEL)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,icpu,i_err)
!#endif
IF(icpu.eq.0)parent=.true.
END
!****************************************************************************************!

SUBROUTINE GlobSumR(myreal_in,myreal_out,leng)
IMPLICIT NONE
!#if defined (_PARALLEL)
INCLUDE 'mpif.h'
!#endif
REAL*8 :: myreal_in(*), myreal_out(*)
INTEGER :: leng,i_err
!#if defined (_PARALLEL)
CALL MPI_Allreduce(myreal_in,myreal_out,leng,MPI_DOUBLE_PRECISION, &
                   MPI_SUM,MPI_COMM_WORLD,i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE GlobSumR_1(myreal_in,myreal_out,leng)
IMPLICIT NONE
!#if defined (_PARALLEL)
INCLUDE 'mpif.h'
!#endif
REAL*8 :: myreal_in, myreal_out
INTEGER :: leng,i_err
!#if defined (_PARALLEL)
CALL MPI_Allreduce(myreal_in,myreal_out,leng,MPI_DOUBLE_PRECISION, &
                   MPI_SUM,MPI_COMM_WORLD,i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE GlobSumI(myint_in,myint_out,leng)
IMPLICIT NONE
!#if defined (_PARALLEL)
INCLUDE 'mpif.h'
!#endif
INTEGER :: myint_in(*), myint_out(*)
INTEGER :: leng,i_err
!#if defined (_PARALLEL)
call MPI_Allreduce(myint_in,myint_out,leng,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,i_err)
!#endif
END
!****************************************************************************************!
! Distributing over MTD Steps
SUBROUTINE Distribute_Mtd_Min_Max(i_md,md_steps,w_cv,w_hill,mtd_min,mtd_max)
IMPLICIT NONE
INTEGER, INTENT(IN)::md_steps, w_cv, w_hill, i_md
INTEGER, INTENT(INOUT):: mtd_min, mtd_max

!Local Variables
INTEGER :: i,ncpu,icpu,mtd_max_old

CALL Get_ncpu(ncpu)
CALL Get_cpuid(icpu)

!if (icpu.eq.0)  write(*,*)'NCPUS',ncpu

! Calulate mtd_max step

! As per shivani code
 mtd_max=((i_md-1)*w_cv/w_hill)+1
 IF( MOD((i_md-1)*w_cv,w_hill).eq.0)  mtd_max=mtd_max-1

! As per Anji Code
!  mtd_max=(i_md*w_cv/w_hill)

!if (icpu.eq.0)print*,"MTD MAX:",mtd_max
! with 1 processor
  if(ncpu.eq.1) then
    mtd_min=1
    mtd_max=mtd_max
  end if

CALL Sync_procs
! with greater than 1 processors (MPI)

IF(ncpu.GT.1) THEN
 mtd_max_old=mtd_max
 mtd_max=mtd_max/ncpu

  CALL Sync_procs

 mtd_min=(mtd_max*icpu)+1

 if(icpu.eq.ncpu-1) then
     mtd_max=mtd_max*(icpu+1)+mod(mtd_max_old,ncpu)
 else
    mtd_max=mtd_max*(icpu+1)
 endif
END IF

END
!****************************************************************************************!

