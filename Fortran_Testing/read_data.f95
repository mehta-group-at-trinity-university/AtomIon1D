program read_data

implicit none

! declaration of variables
real :: x,y,z

!main

open(10,file='sample_data.txt')

read(10,*) x,y,z

print *,x,y,z


end program read_data
