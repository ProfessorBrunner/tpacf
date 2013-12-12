diff -q atest.dat_ddbins answers/add_answer > alog 2>> alog
diff -q atest.dat_atest.dat_drbins answers/adr_answer >> alog 2>> alog
diff -q atest.dat_rrbins answers/arr_answer >> alog 2>> alog
echo "Problems detected in angular calculation:"
wc -l < alog
diff -q stest.dat_ddbins answers/sdd_answer > slog 2>> slog
diff -q stest.dat_stest.dat_drbins answers/sdr_answer >> slog 2>> slog
diff -q stest.dat_rrbins answers/srr_answer >> slog 2>> slog
echo "Problems detected in spatial calculation:"
wc -l < slog
