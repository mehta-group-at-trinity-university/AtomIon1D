CMP     = gfortran
CMPFLAGS = -ffixed-line-length-132 -O3
DEBUG   = -fcheck=all
FORCEDP = #-fdefault-real-8 -fdefault-double-8
INCLUDE =  -I/opt/homebrew/include/
LAPACK =  -framework accelerate
ARPACK =  -I/usr/local/include -L /usr/local/lib/ -larpack
OBJS  = besselnew.o Bsplines.o matrix_stuff.o RMATPROP2016.o Quadrature.o

RMATPROP2016.x:	   ${OBJS}
	${CMP} ${DEBUG} ${OBJS} ${INCLUDE} ${ARPACK} ${LAPACK}  ${CMPFLAGS} ${FORCEDP} -o RMATPROP2016.x

RMATPROP2016.o: RMATPROP2016.f90
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c RMATPROP2016.f90

matrix_stuff.o: matrix_stuff.f
	${CMP} ${FORCEDP} ${CMPFLAGS} -c matrix_stuff.f

Bsplines.o:	Bsplines.f
	${CMP} ${FORCEDP} ${CMPFLAGS} -c Bsplines.f

nrtype.mod: modules_qd.o
	${CMP} ${FORCEDP} modules_qd.o

modules_qd.o:	modules_qd.f90
	${CMP} ${FORCEDP} -c modules_qd.f90

besselnew.o:	besselnew.f
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c besselnew.f

Quadrature.mod: Quadrature.o
		${CMP} ${FORCEDP} Quadrature.o

Quadrature.o: Quadrature.f90
	${CMP} ${FORCEDP} -c Quadrature.f90

clean:
	rm -f *.mod *.o *.x


HHL1DHyperspherical.x: modules_qd.o besselnew.o Bsplines.o  ../bspllib_22/bspline90_22.o
	gfortran -O -ffixed-line-length-132 -fcheck=all -Wmaybe-uninitialized Bsplines.o modules_qd.o bspline90_22.o besselnew.o -L/opt/homebrew/lib/ -larpack -framework accelerate  HHL1DHyperspherical.Working.f -o HHL1DHyperspherical.x


Bsplines.o:	Bsplines.f
	gfortran -O  -ffixed-line-length-132 -c Bsplines.f	

besselnew.o :	besselnew.f
	gfortran -O  -ffixed-line-length-132 -c besselnew.f	

modules_qd.o :	modules_qd.f90
	gfortran -O  -c modules_qd.f90

bspline90_22.o:	 bspline90_22.f90
	gfortran -O -c bspline90_22.f90

