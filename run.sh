#!/bin/bash

/opt/intel/impi/4.0.3.008/intel64/bin/mpicc /scratch/apc/std202013/$2 -o output
/opt/intel/impi/4.0.3.008/intel64/bin/mpirun -np $1 /scratch/apc/std202013/output

