program helloworld

  double precision x,y,z
  integer f,fold
  write(6,*) "Hello World!"
  write(6,*) "My name is Rosa!"


  x = 20d0
  y = 3d0

  z = x*y
  write(6,*) "z = ", z

  f = 0
  do n = 1,20
     f = fold+n
     write(6,*) f
     fold = f
     
  enddo
     
end program helloworld
