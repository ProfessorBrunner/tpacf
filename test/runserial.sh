./cleartest.sh
cp ../bin/* .
cp testdata/* .
echo
./precompute afilelist 0 8
echo
./precompute sfilelist 1 8
echo
echo "Computing angular correlation"
./correlate aparams.in
echo
echo "Computing spatial correlation"
./correlate sparams.in
echo
rm *.dat *filelist *params.in
./compare.sh
