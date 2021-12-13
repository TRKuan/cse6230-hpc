# We want as reproducible an environment as possible
module purge
module load git nvhpc

# clean temporary files
make clean

make main

echo "["

FIRST_ITER="yes"
for inputs in "100000 5" "1000 10" "100 15" "10 20" "10 25" "10 30" "10 35" "10 40" "10 45" "10 50" "10 55" "10 60"; do
  if [ "$FIRST_ITER" = "yes" ]; then
    FIRST_ITER="no"
  else
    echo ","
  fi
  ./main float $inputs
  echo ","
  ./main double $inputs
done

echo "]"


