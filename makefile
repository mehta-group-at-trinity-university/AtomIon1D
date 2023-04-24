CMP     = gfortran
CMPFLAGS = -ffixed-line-length-132 -O3
DEBUG   = -fcheck=all
FORCEDP = #-fdefault-real-8 -fdefault-double-8
INCLUDE =  -I/usr/local/
LAPACK =  -framework Accelerate
ARPACK =  -I/usr/local/include -L /usr/local/lib/ -larpack
OBJS  = besselnew.o Bsplines.o Quadrature.o AtomIon1D.o
#matrix_stuff.o

AtomIon1D.x:	   ${OBJS}
	${CMP} ${DEBUG} ${OBJS} ${INCLUDE} ${ARPACK} ${LAPACK}  ${CMPFLAGS} ${FORCEDP} -o AtomIon1D.x

AtomIon1D.o: AtomIon1D.f
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c AtomIon1D.f

#matrix_stuff.o: matrix_stuff.f
#	${CMP} ${FORCEDP} ${CMPFLAGS} -c matrix_stuff.f

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

