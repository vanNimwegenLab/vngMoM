#!/bin/bash
#SBATCH --qos=6hours
#SBATCH --mem-per-cpu=2G
#SBATCH -o slogs/$JOB_NAME.o$JOB_ID
#SBATCH -e slogs/$JOB_NAME.e$JOB_ID

# syntax: get_size_and_fluo.sh ./path/to/get_size_and_fluo_basic.pl ./path/to/input ./path/to/output

INPUT=$(readlink -f $2)
INPUT_FILENAME=$(basename "$INPUT")
INPUT_DIR=$(dirname "$INPUT")
OUTPUT=$(readlink -m $3) # -m required because the path probably doesn't exist yet
DIR=$TMPDIR

# move to a tmp location and copy the input file (erik's scripts work only in the working dir)
cd $DIR
cp $INPUT .

# run perl script and extract output file name
OUT=$(perl $1 $INPUT_FILENAME)
OUTNAME=$(echo $OUT | perl -pe 's/^>//') # replace `^>` with `` (using regex)

# housekeeping
mkdir -p $(dirname "$OUTPUT") # -p: no error if existing, make parent directories as needed
mv $OUTNAME $OUTPUT
rm $INPUT_FILENAME
