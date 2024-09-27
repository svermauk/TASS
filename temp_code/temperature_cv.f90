program temperature_cv

implicit none

real*8, parameter :: gasconst = 8.314E-3 ! in kJ mol^-1 K^-1
real*8, parameter :: pi=4.d0*atan(1.d0)
real*8, parameter :: convert_to_mass=100.d0   ! changing tau to mass conversion where mass is in amu A^-2 d^-2
real*8, parameter :: convert_to_K=0.01
real*8, allocatable :: tau(:), kappa(:), mass(:)
real*8 ::temp, temp_av
integer :: ncv, istep, ios,i
character :: filename*20

real*8, allocatable :: vfict(:)
open(1,file='input',status='old')
write(*,*)'number of CVs'
read(1,*)ncv

allocate(tau(ncv))
allocate(mass(ncv))
allocate(kappa(ncv))

do i= 1,ncv
write(*,*) 'reading details about CV number =', i
!write(*,*)'fetching the value of tau (as in plumed)'
read(1,*) tau(i)

!write(*,*)'fetching the value of kappa (as in plumed)'
read(1,*) kappa(i)

mass(i)= tau(i)*tau(i)*kappa(i)/4.d0/pi/pi*convert_to_mass
print *, "mass in amu Ang.^-2  d^-2 : ", mass(i)

end do

!convert tau to mass units (i.e. a.m.u. Angstrom^-2 d^-2) ; here d is the units of CV


allocate(vfict(ncv))

write(*,*)'enter the plumed ouput file name with vfict'
read(*,*)filename

open(1,file=trim(filename),status='old',form='formatted')
open(2,file='cv_temp_test.dat')
print *, 'NOTE: this program will skip the first line of the vfict file'
print *, 'NOTE: this program skips the first column of the vfict file'
read(1,*)

istep=0
temp_av=0.d0
do 
  read(1,*,IOSTAT=ios)temp,vfict(1:ncv) !temp is dummy here; vfict is the velocity in d/ps units
  if(ios.ne.0)exit
  istep=istep+1
  temp=0.d0
  do i=1,ncv
!  temp=temp+mass(i)*dot_product(vfict(i),vfict(i))/gasconst*convert_to_K/dfloat(ncv)
  temp=temp+mass(i)*vfict(i)*vfict(i)/gasconst*convert_to_K
  end do
  temp = temp/dfloat(ncv)
  temp_av=temp_av+temp
  print *, istep,temp,temp_av/dfloat(istep)
  write(2,*) istep,temp,temp_av/dfloat(istep)
end do
print *, 'total md steps read from output =',istep
print *, 'output writtten in cv_temp_test.dat'
close(1)
close(2)

end program temperature_cv
