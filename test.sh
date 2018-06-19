#############################################################################
#
#  test.sh              Andrés Herrera-Poyatos <andreshp9@gmail.com>
#                       Pedro A. Garcia-Sanchez <pedro@ugr.es>
#
# Script which launch the tests for the kronecker algorithms.
#
#############################################################################

mkdir -p data
mkdir -p images

#echo "Executing the first test for all the algorithms..."
#bash ./tests/test_kronecker_all_1.sh | tail -n +5  > ./data/all_1.csv
#echo "Executing the second test for all the algorithms..."
#bash ./tests/test_kronecker_all_2.sh | tail -n +5  > ./data/all_2.csv
#echo "Executing the third test for all the algorithms..."
#bash ./tests/test_kronecker_all_pedro.sh | tail -n +5  > ./data/all_pedro.csv
#echo "Executing the first test for the fastest algorithms..."
#bash ./tests/test_kronecker_fast_1.sh | tail -n +5  > ./data/fast_1.csv
#echo "Executing the second test for the fastest algorithms..."
#bash ./tests/test_kronecker_fast_2.sh | tail -n +5  > ./data/fast_2.csv
#echo "Executing the third test for the fastest algorithms..."
#bash ./tests/test_kronecker_fast_3.sh | tail -n +5  > ./data/fast_3.csv
echo "Executing the fourh test for the fastest algorithms..."
bash ./tests/test_kronecker_fast_4.sh | tail -n +5  > ./data/fast_4.csv
#echo "Executing the fith test for all the algorithms..."
#bash ./tests/test_kronecker_fast_pedro.sh | tail -n +5  > ./data/fast_pedro.csv
#pl ./data/all_1.csv -t "Kronecker algorithms on x^{2^n}-1" -x "n" -y "Time ms" -l -o images/all_1.png
#pl ./data/all_2.csv -t "Kronecker algorithms on x^{2^n}+2" -x "n" -y "Time ms" -l -o images/all_2.png
#pl ./data/all_pedro.csv -t "Kronecker algorithms on 1-x+x^{n}-x^{2n-1}+x^{2n}" -x "n" -y "Time ms" -l -o images/all_pedro.png
pl ./data/fast_1.csv -t "Kronecker algorithms on x^{2^n}-1" -x "n" -y "Time ms" -l -o images/fast_1.png
pl ./data/fast_2.csv -t "Kronecker algorithms on x^{2 2^n}+x^{2^n}-2" -x "n" -y "Time ms" -l -o images/fast_2.png
pl ./data/fast_3.csv -t "Kronecker algorithms on Phi_{2^n}" -x "n" -y "Time ms" -l -o images/fast_3.png
pl ./data/fast_4.csv -t "Kronecker algorithms on x^{2^n}+2" -x "n" -y "Time ms" -l -o images/fast_4.png
pl ./data/fast_pedro.csv -t "Kronecker algorithms on 1-x+x^{n}-x^{2n-1}+x^{2n}" -x "n" -y "Time ms" -l -o images/fast_pedro.png

