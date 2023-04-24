1D+.x:	bspline_22.o 1D+.o
	gfortran	 bspline_22.o 1D+.o  -framework Accelerate -L/usr/local/lib/ -larpack

check1:
	echo "check 1"

1D+.o:	1D+.f
	gfortran	-ffixed-line-length-132  -c 1D+.f

bspline_22.o:	bspline_22.f
	gfortran	-ffixed-line-length-132 -c bspline_22.f

clean:
	rm -f *.o *.x fort.*

checkf:
	@echo "make file deleted files probably"