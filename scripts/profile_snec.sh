# Run snec with varying number of cores and output to file timing

#!/bin/bash

Ncores=(1 5 10 15 20 30)

for i in "${Ncores[@]}"; do
    export OMP_NUM_THREADS=$i
    {
	echo $OMP_NUM_THREADS >&2
	time $SNEC_DIR/snec
    } 2>> $SNEC_DIR/output_profiling.txt
done
