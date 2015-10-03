#!/bin/csh
unset noclobber

#--- initialize test mode
touch RSTTEST

#--- remove files from possible previous test
rm -f quad_rstini
rm -f one*
rm -f two*
rm -f full*

#--- set namelist parameter nsteps in quad_namelist to 10000 
sed -i -e 's/nsteps.*/nsteps = 10000/' quad_namelist

#--- compile and run first part of run
make
make run
cp quad_diag   one_quad_diag
cp quad_gp     one_quad_gp
cp quad_rstfin one_quad_rstfin
cp quad_tseri  one_quad_tseri
mv quad_rstfin quad_rstini
#--- run second part of run
make run
cp quad_diag   two_quad_diag
cp quad_gp     two_quad_gp
cp quad_rstfin two_quad_rstfin
cp quad_tseri  two_quad_tseri

#--- prepare namelist file for total run
sed -i -e 's/nsteps.*/nsteps = 20000/' quad_namelist

#--- prepare and do total run
rm quad_rstini
make run
cp quad_diag   full_quad_diag
cp quad_gp     full_quad_gp
cp quad_rstfin full_quad_rstfin
cp quad_tseri  full_quad_tseri

#--- compare if final restart files differ
diff two_quad_rstfin full_quad_rstfin

if ($? == 0) echo "Restart test passed"
