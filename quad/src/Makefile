NGX=64

F90=gfortran
OPTS=-O2 -ffpe-trap=invalid,zero,overflow  -finit-real=snan
DEBOPTS=-g -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -fcheck=all -finit-real=snan

%.o : %.c
	$(CC) -c $<

%.o : %.f90
	$(F90) -c $(OPTS) $<

quad.x:	quad.o fftfor.o fftasm.o rstmod.o
	$(F90) -o quad.x $(OPTS) quad.o fftfor.o fftasm.o rstmod.o
quad.o:	quad.f90
fftfor.o:	fftfor.f90

fftasm.o:	fftasm.s
	as -o fftasm.o fftasm.s

rstmod.o:	rstmod.f90

clean:
	rm -f quad.x *.o *.mod

cleanall:
	rm -f quad.x *.o *.mod quad_gp quad_gp.nc quad_tseri quad_cfl quad_diag quad_rstini quad_rstfin

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
allrst:
	make cleanall
	make
	make run
	make nc
	cp quad_gp.nc quad_gp_01.nc
	cp quad_rstfin quad_rstini
	make run
	make nc
	cp quad_gp.nc quad_gp_02.nc
