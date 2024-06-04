CMP     = gfortran
CMPFLAGS = -ffixed-line-length-132 -O3
DEBUG   = -fcheck=all
FORCEDP = #-fdefault-real-8 -fdefault-double-8
ARPACK =  -L/opt/homebrew/lib/ -larpack
INCLUDE =  -I/opt/homebrew/include
LAPACK =  -framework Accelerate
OBJS  = besselnew.o Bsplines.o Quadrature.o noninteracting1D.o
#matrix_stuff.o

noninteracting1D.x:	   ${OBJS}
	${CMP} ${DEBUG} ${OBJS} ${INCLUDE} ${ARPACK} ${LAPACK}  ${CMPFLAGS} ${FORCEDP} -o noninteracting1D.x

noninteracting1D.o: noninteracting1D.f90
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c noninteracting1D.f90

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

