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


echo "Executing the first test for all the algorithms..."
bash ./tests/test_kronecker_all_1.sh | tail -n +5  > ./data/all_1.csv
echo "Executing the second test for all the algorithms..."
bash ./tests/test_kronecker_all_2.sh | tail -n +5  > ./data/all_2.csv
pl ./data/all_1.csv -t "Kronecker algorithms on x^{2^n}-1" -x "n" -y "Time ms" -l -o images/all_1.png
pl ./data/all_2.csv -t "Kronecker algorithms on x^{2^n}+2" -x "n" -y "Time ms" -l -o images/all_2.png

