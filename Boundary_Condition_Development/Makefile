1D+.x:	bspline_22.o 1D+.o
	gfortran	 bspline_22.o 1D+.o  -framework Accelerate -L/usr/local/lib/ -larpack

1D+.o:	1D+.f
	gfortran	-ffixed-line-length-132  -c 1D+.f

check1:
	@echo "1 make file deleted files"

bspline_22.o:	bspline_22.f
	gfortran	-ffixed-line-length-132 -c bspline_22.f

check2:
	@echo "2 make file deleted files"

clean:
	rm -f *.o *.x fort.*

check3:
	@echo "3 make file deleted files"