1D+.x:	bspline_22.o 1D+.o
	gfortran	 bspline_22.o 1D+.o  -framework accelerate -L/usr/local/lib/ -larpack

1D+.o:	1D+.f
	gfortran	-ffixed-line-length-132  -c 1D+.f

bspline_22.o:	bspline_22.f
	gfortran	-ffixed-line-length-132 -c bspline_22.f

clean:
	rm -f *.o *.x fort.*
