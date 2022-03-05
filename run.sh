#$ -cwd
#$ -S /bin/bash
#$ -N run
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID

# Hrothgar Cluster
# $ -q west
# $ -pe west 12
# $ -P hrothgar

# Quanah Cluster
#$ -q omni
#$ -pe fill 3
#$ -P quanah

# mpirun --machinefile machinefile.$JOB_ID -np $NSLOTS 
 ./source 2 2250 80
# ./source -dim 33 -ndata 68040  -file ColorHistogram.dat
# ./source -dim 10 -ndata 45730  -file CASP.csv.bin
# ./source -dim 29 -ndata 11000000  -file HIGGS.dat


