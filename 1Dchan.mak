1Dchan.x:	1Dchan.f Bsplines.f matrix_stuff.o modules_qd.o besselnew.o 1Dchan.o Bsplines.o  1Dpot.o ../bspllib_22/bspline90_22.o
	gfortran -O  1Dchan.o  Bsplines.o 1Dpot.o modules_qd.o ../bspllib_22/bspline90_22.o besselnew.o matrix_stuff.o  -L/usr/local/lib/ -larpack -framework Accelerate  -lm -o 1Dchan.x

1Dpot.o:	1Dpot.f
	gfortran -O  -ffixed-line-length-132 -c 1Dpot.f

1Dchan.o:	1Dchan.f
	gfortran -O  -ffixed-line-length-132 -c 1Dchan.f

Bsplines.o:	Bsplines.f
	gfortran -O  -ffixed-line-length-132 -c Bsplines.f	

besselnew.o :	besselnew.f
	gfortran -O  -ffixed-line-length-132 -c besselnew.f	

modules_qd.o :	modules_qd.f90
	gfortran -O  -c modules_qd.f90	

matrix_stuff.o:	matrix_stuff.f
	gfortran -O -ffixed-line-length-132 -c matrix_stuff.f
