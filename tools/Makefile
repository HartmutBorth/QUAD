NGX=64

F90=gfortran

%.o : %.c
	$(CC) -c $<

%.o : %.f90
	$(F90) -c $(OPTS) $<

MakeIni.x:	MakeIni.o
	$(F90) -o MakeIni.x MakeIni.f90

clean:
	rm -f MakeIni.x *.o *.mod

cleanall:
	rm -f MakeIni.x *.o *.mod G*.srv

run:
	MakeIni.x $(NGX)

all:
	make cleanall
	make
	make run
