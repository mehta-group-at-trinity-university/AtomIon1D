program write_data

implicit none

! declaration of variables
real :: x,y,z

x = 1
y = 2
z = 3

!main

open(12,file='written_data.txt')

write(12,*) x,y,z

print *, 'Nice data'

end program write_data
