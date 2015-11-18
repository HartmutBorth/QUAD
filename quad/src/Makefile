NGX=64

F90=gfortran
OPTS=-O2 -ffpe-trap=invalid,zero,overflow  -finit-real=snan
DEBOPTS=-g -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -fcheck=all -finit-real=snan

%.o : %.c
	$(CC) -O3 -c $<

%.o : %.f90
	$(F90) -c $(OPTS) $<

quad.x:	quad.o fftfor.o fftasm.o rstmod.o guimod.o pumax.o
	$(F90) -o quad.x $(OPTS) quad.o fftfor.o fftasm.o rstmod.o guimod.o pumax.o -L/opt/X11 -lX11
quad.o:	quad.f90
fftfor.o:	fftfor.f90

fftasm.o:	fftasm.s
	as -o fftasm.o fftasm.s

rstmod.o:	rstmod.f90

guimod.o:	guimod.f90

pumax.o:	pumax.c

clean:
	rm -f quad.x *.o *.mod

cleanall:
	rm -f quad.x *.o *.mod quad_gp quad_gp.nc quad_sp quad_sp.nc quad_tseri quad_cfl quad_diag quad_rstini quad_rstfin
	make cleanrsttest

run:
	time quad.x $(NGX)

nc:	
	cdo -f nc copy quad_gp quad_gp.nc

view:
	ncview quad_gp.nc
all:
	make cleanall
	make
	make run
	make nc
	make view
rsttest:
	cp ../../tools/quad_rst_test.csh .       
	csh quad_rst_test.csh
	rm -f quad_rst_test.csh
cleanrsttest:
	cp ../../tools/clean_quad_rst_test.csh .       
	csh clean_quad_rst_test.csh
	rm -f clean_quad_rst_test.csh
