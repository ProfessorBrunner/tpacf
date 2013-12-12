cp ../bin/* .
cp testdata/* .
./aprecompute afilelist 0 8
./sprecompute sfilelist 1 8
bsub < cmpi_angular.lsf
bsub < cmpi_spatial.lsf
rm *.dat *filelist 
