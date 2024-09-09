#!/bin/bash
#SBATCH --job-name=qc_plots
#SBATCH --output=test.out
#SBATCH --error=test.err

touch test_file
echo "Testing directory access" > test_output.txt
