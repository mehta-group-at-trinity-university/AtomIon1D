#!/bin/bash

# IN THIS SCRIPT, A WORKING DIRECTORY IS CREATED
# SO THAT MULTIPLE RUNS CAN HAPPEN IN PARALLEL.

# NOTE THAT NewFitProg.f and Legendre.dat SHOULD
# BE MOVED TO CABA, WHERE NewFitProg.f IS COMPILED.
# THEN, NewFitProg.x IS COPIED TO THE WORKING DIR
# (AND PROBABLY ALSO Legendre.dat). 
# THIS WAY, IT SHOULD BE CLEANER


BOUNDSTATES=1
PARITY=P
OPENCHANNELS=1
# these parameter will appear in directories and file names


AIDIR=${BOUNDSTATES}bs_R200
# subdirectory where AtomIon outputs are stored

PARFILE=CABA_${PARITY}${OPENCHANNELS}.inp
# name of input file


# —————————————————————————————————————————————————————
# things below this line should not need to be modified

NORMAL=$(tput sgr0)
RED=$(tput setaf 1)
BLUE=$(tput setaf 4)
GREEN=$(tput setaf 2)
# colors for printed text


CABADIR=${BOUNDSTATES}bs 
# subdirectory where CABA files are stored

WORKDIR=`mktemp -d -p "$CABADIR"`
# temporary working directory with random name


if [ ! -d "$WORKDIR" ]
then
  echo "${RED}Could not create working directory.${NORMAL}"
  exit 1
else
  echo "Working directory ${BLUE}$WORKDIR${NORMAL} created;"
  sleep 1
fi
# check if directory was created


cp ../../data_local/$AIDIR/DiabaticEnergies.dat ./$WORKDIR
cp ../../data_local/$AIDIR/AdiabaticEnergies.dat ./$WORKDIR
cp ../../data_local/$AIDIR/VQmat.dat ./$WORKDIR
cp ../../data_local/$AIDIR/Pmat.dat ./$WORKDIR
cp ../../data_local/$AIDIR/QuickEffV.dat ./$WORKDIR
cp ../../data_local/$AIDIR/QuickPMat.dat ./$WORKDIR
# copies output of AtomIon1D from data_local 
# to working subdirectory in CABA directory 

if [ $? -ne 0 ]
then
    echo "${RED}No AtomIon data for these parameters.${NORMAL}"
    exit 1
fi

make
# compiles CABA code

if [ $? -ne 0 ]
then
    echo "${RED}MyCABA does not compile.${NORMAL}"
    exit 1
fi

cp MyCABA.EXT.x ./$WORKDIR
echo "CABA executable copied to working directory in ${BLUE}$CABADIR${NORMAL}"
# copies executable to bound-state subdirectory

gfortran -ffixed-line-length-none -O3 NewFitProg.f -o NewFitProg.x
# compiles fit

cp NewFitProg.x ./$WORKDIR
cp Legendre.dat ./$WORKDIR
echo "Fit executable and Legendre file copied to working directory in ${BLUE}$CABADIR${NORMAL}"
# copies fit executable and Legendre file to working subdirectory

cd ./$WORKDIR
pwd
# goes to bound-state subdirectory

sleep 1

if [ -f comment.txt ]
then 
  rm comment.txt
fi

echo "Do you want to leave a comment?"
select yn in "Yes" "No"; 
do 
  case $yn in
      Yes ) echo "Write comment and press enter:";
            read comment
            echo "$comment" >> comment.txt
            break;;
      No  ) echo "No comment." >> comment.txt
            break;;
  esac
done
# gives you the possibility to write a comment in a txt file

./NewFitProg.x
# runs fit program

if [ $? -ne 0 ]
then
    echo "${RED}Problem with fit.${NORMAL}"
    exit 1
fi

./MyCABA.EXT.x < ../input_files/$PARFILE
# runs CABA with the specified input file

if [ $? -ne 0 ]
then
    echo "${RED}Some problem occurred while running MyCABA.${NORMAL}"
    exit 1
fi

echo "${GREEN}CABA code complete;${NORMAL}"

sleep 1

rm AdiabaticEnergies.dat 
rm DiabaticEnergies.dat
rm Pmat.dat 
rm VQmat.dat 
rm QuickEffV.dat 
rm QuickPMat.dat  
# removes AtomIon1D output from working directory

check=0
i=1
while [ $check == 0 ]
do
    DIRECTORY=../../../../data_local/$AIDIR/${PARITY}${OPENCHANNELS}_set${i}
    if [ -d "$DIRECTORY" ]
    then
        i=$((i+1))
    else
        check=1
    fi
done

DATADIR=${PARITY}${OPENCHANNELS}_set${i}

mkdir ../../../../data_local/$AIDIR/$DATADIR
# creates new subdirectory to avoid overwriting data

cp ../input_files/$PARFILE ../../../../data_local/$AIDIR/$DATADIR/
# copies parameter file in data subdirectory

mv QuickVQMat.dat ../../../../data_local/$AIDIR/
mv comment.txt ../../../../data_local/$AIDIR/$DATADIR/
mv Kmatrix*.dat ../../../../data_local/$AIDIR/$DATADIR/
mv Observables*.dat ../../../../data_local/$AIDIR/$DATADIR/
mv Smatrix*.dat ../../../../data_local/$AIDIR/$DATADIR/
mv TimeDelay*.dat ../../../../data_local/$AIDIR/$DATADIR/
mv grid.used.data ../../../../data_local/$AIDIR/$DATADIR/
# moves all output in data subdirectory

if [ $? -ne 0 ]
then
    echo "${RED}Some problem occurred while moving the files.${NORMAL}"
    exit 1
else 
    echo "Output moved subdirectory ${BLUE}$AIDIR${NORMAL} in ${BLUE}data_local${NORMAL};"
    cd ../..
    rm -rf $WORKDIR
    echo "Working directory deleted;"
    echo "Number of open channels is ${BLUE}$OPENCHANNELS${NORMAL}, set number is ${BLUE}$i${NORMAL}." 
fi

