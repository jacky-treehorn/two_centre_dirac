FORTRAN=gfortran
FLAG=-O2
FLAG_TST=-check all -warn unused
COMPILE=-c
OUTPUT=-o
ROLLING=-funroll-loops
TRACING=-traceback
DEBUGGING=-g

$(shell mkdir -p objects)
$(shell mkdir -p binaries)

binaries/B_spline.exe:	objects/new_splines.o objects/common_subroutines.o objects/main.o objects/potential.o objects/io.o objects/my_blas.o
	$(FORTRAN) $(FLAG) $(ROLLING) $(OUTPUT) B_spline.exe objects/new_splines.o objects/common_subroutines.o objects/main.o objects/potential.o objects/io.o objects/my_blas.o
	mv B_spline.exe binaries
objects/new_splines.o:	sources/inc.par sources/new_splines.f
	$(FORTRAN) $(FLAG) $(ROLLING) $(COMPILE) sources/new_splines.f
	\mv *.o objects
objects/common_subroutines.o:	sources/inc.par sources/common_subroutines.f
	$(FORTRAN) $(FLAG) $(ROLLING) $(COMPILE) sources/common_subroutines.f
	\mv *.o objects
objects/main.o:	sources/inc.par sources/main.f
	$(FORTRAN) $(FLAG) $(ROLLING) $(COMPILE) sources/main.f
	\mv *.o objects
objects/potential.o:	sources/inc.par sources/potential.f
	$(FORTRAN) $(FLAG) $(ROLLING) $(COMPILE) sources/potential.f
	\mv *.o objects
objects/io.o:	sources/io.f
	$(FORTRAN) $(FLAG) $(ROLLING) $(COMPILE) sources/io.f
	\mv *.o objects
objects/my_blas.o:	sources/my_blas.f
	$(FORTRAN) $(FLAG) $(ROLLING) $(COMPILE) sources/my_blas.f
	\mv *.o objects
binaries/B_spline_tst.exe:	objects/new_splines_tst.o objects/common_subroutines_tst.o objects/main_tst.o objects/potential_tst.o objects/io_tst.o objects/my_blas_tst.o
	$(FORTRAN) $(FLAG_TST) $(OUTPUT) B_spline_tst.exe objects/new_splines_tst.o objects/common_subroutines_tst.o objects/main_tst.o objects/potential_tst.o objects/io_tst.o objects/my_blas_tst.o
	mv B_spline_tst.exe binaries
objects/new_splines_tst.o:	sources/inc.par sources/new_splines.f
	$(FORTRAN) $(FLAG_TST) $(COMPILE) sources/new_splines.f
	\mv new_splines.o objects/new_splines_tst.o
objects/common_subroutines_tst.o:	sources/inc.par sources/common_subroutines.f
	$(FORTRAN) $(FLAG_TST) $(COMPILE) sources/common_subroutines.f
	\mv common_subroutines.o objects/common_subroutines_tst.o
objects/main_tst.o:	sources/inc.par sources/main.f
	$(FORTRAN) $(FLAG_TST) $(COMPILE) sources/main.f
	\mv main.o objects/main_tst.o
objects/potential_tst.o:	sources/inc.par sources/potential.f
	$(FORTRAN) $(FLAG_TST) $(COMPILE) sources/potential.f
	\mv potential.o objects/potential_tst.o
objects/io_tst.o:	sources/io.f
	$(FORTRAN) $(FLAG_TST) $(COMPILE) sources/io.f
	\mv io.o objects/io_tst.o
objects/my_blas_tst.o:	sources/my_blas.f
	$(FORTRAN) $(FLAG_TST) $(COMPILE) sources/my_blas.f
	\mv my_blas.o objects/my_blas_tst.o
binaries/B_spline_g.exe:	objects/new_splines_g.o objects/common_subroutines_g.o objects/main_g.o objects/potential_g.o objects/io_g.o objects/my_blas_g.o
	$(FORTRAN) $(DEBUGGING) $(OUTPUT) B_spline_g.exe objects/new_splines_g.o objects/common_subroutines_g.o objects/main_g.o objects/potential_g.o objects/io_g.o objects/my_blas_g.o
	mv B_spline_g.exe binaries
objects/new_splines_g.o:	sources/inc.par sources/new_splines.f
	$(FORTRAN) $(DEBUGGING) $(COMPILE) sources/new_splines.f
	\mv new_splines.o objects/new_splines_g.o
objects/common_subroutines_g.o:	sources/inc.par sources/common_subroutines.f
	$(FORTRAN) $(DEBUGGING) $(COMPILE) sources/common_subroutines.f
	\mv common_subroutines.o objects/common_subroutines_g.o
objects/main_g.o:	sources/inc.par sources/main.f
	$(FORTRAN) $(DEBUGGING) $(COMPILE) sources/main.f
	\mv main.o objects/main_g.o
objects/potential_g.o:	sources/inc.par sources/potential.f
	$(FORTRAN) $(DEBUGGING) $(COMPILE) sources/potential.f
	\mv potential.o objects/potential_g.o
objects/io_g.o:	sources/io.f
	$(FORTRAN) $(DEBUGGING) $(COMPILE) sources/io.f
	\mv io.o objects/io_g.o
objects/my_blas_g.o:	sources/my_blas.f
	$(FORTRAN) $(DEBUGGING) $(COMPILE) sources/my_blas.f
	\mv my_blas.o objects/my_blas_g.o
pu:	
	\rm *~ */*~ */*/*~ */*/*/*~
