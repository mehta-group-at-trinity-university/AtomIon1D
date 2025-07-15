CMP = gfortran
STD = 
CMPFLAGS = -ffixed-line-length-none -O3
DEBUG   = -fcheck=all
FORCEDP = #-fdefault-real-8 -fdefault-double-8
ARPACK =  -L/opt/homebrew/lib/ -larpack
INCLUDE =  -I/opt/homebrew/include
LAPACK =  -framework Accelerate
OBJS  =  akima.o eispack.o OurCABA.EXT.o

OurCABA.EXT.x:  ${OBJS}
	${CMP} ${DEBUG} ${OBJS} ${INCLUDE} ${ARPACK} ${LAPACK}  ${CMPFLAGS} ${FORCEDP} -o OurCABA.EXT.x

OurCABA.EXT.o: OurCABA.EXT.f
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c OurCABA.EXT.f

eispack.o: eispack.f
	${CMP} ${STD} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c eispack.f

akima.o: akima.f
	${CMP} ${STD} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c akima.f
clean:
	rm -f *.mod *.o *.x
