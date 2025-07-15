#!/bin/bash

TOT_BOUNDSTATES=3
# number of total boundstates considered in the AtomIon code

USED_BOUNDSTATES=2  # number of bound states to consider
COLLISIONAL=12      # number of collisional channels
PARITY=P
OPENCHANNELS=1
# these parameter will appear in directories and file names

AIDIR=${TOT_BOUNDSTATES}bs_1200xPoints_${COLLISIONAL}collisional
# subdirectory where AtomIon outputs are stored

TOT_STATES=$(($USED_BOUNDSTATES+$USED_BOUNDSTATES+$COLLISIONAL))

PARFILE=OurCABA_${PARITY}${OPENCHANNELS}_${TOT_STATES}channels.inp
# name of input file


# —————————————————————————————————————————————————————
# things below this line should not need to be modified

NORMAL=$(tput sgr0)
RED=$(tput setaf 1)
BLUE=$(tput setaf 4)
GREEN=$(tput setaf 2)
# colors for printed text


CABADIR=${TOT_BOUNDSTATES}bs 
# subdirectory where CABA files are stored

# WORKDIR=`mktemp -d -p "$CABADIR"`
# temporary working directory with random name
# note: with this sintax it is ./CABADIR/WORKDIR


if [ ! -d "../../data_local/$AIDIR" ] 
then
    echo "${RED}Directory${NORMAL} ${BLUE}$AIDIR${NORMAL} ${RED}does not exist.${NORMAL}"
    exit 1
fi

if [ ! -f "./$CABADIR/input_files/$PARFILE" ] 
then
    echo "${RED}Parameter file${NORMAL} ${BLUE}$PARFILE${NORMAL} ${RED}does not exist.${NORMAL}"
    exit 1
fi
# checks if the names of atom-ion directory and CABA parameter file make sense


check=0
i=1
while [ $check == 0 ]
do
    CHECKDIR=./${CABADIR}/workdir${i}
    if [ -d "$CHECKDIR" ]
    then
        i=$((i+1))
    else
        check=1
    fi
done

WORKDIR=workdir${i}
mkdir ./$CABADIR/$WORKDIR
echo "Working directory ${BLUE}$i${NORMAL} created;"


cp ../../data_local/$AIDIR/DiabaticEnergies.dat ./$CABADIR/$WORKDIR
cp ../../data_local/$AIDIR/AdiabaticEnergies.dat ./$CABADIR/$WORKDIR
cp ../../data_local/$AIDIR/VQmat.dat ./$CABADIR/$WORKDIR
cp ../../data_local/$AIDIR/Pmat.dat ./$CABADIR/$WORKDIR
# cp ../../data_local/$AIDIR/QuickEffV.dat ./$CABADIR/$WORKDIR
# cp ../../data_local/$AIDIR/QuickPMat.dat ./$CABADIR/$WORKDIR
cp ../../data_local/$AIDIR/ParAndThresh.dat ./$CABADIR/$WORKDIR
# copies output of AtomIon1D from data_local 
# to working subdirectory in CABA directory 

cp ./$CABADIR/input_files/$PARFILE ./$CABADIR/$WORKDIR

cp Legendre.dat ./$CABADIR/$WORKDIR


make
# compiles CABA code

if [ $? -ne 0 ]
then
    echo "${RED}OurCABA does not compile.${NORMAL}"
    exit 1
fi


cp OurCABA.EXT.x ./$CABADIR/$WORKDIR
echo "CABA executable copied to working directory in ${BLUE}$CABADIR${NORMAL}"
# copies executable to bound-state subdirectory

cd ./$CABADIR/$WORKDIR
pwd
# goes to bound-state subdirectory

sleep 1

if [ -f comment.txt ]
then 
  rm comment.txt
fi

echo "Do you want to leave a comment? [Type number, then press ENTER]"
select yn in "Yes" "No"; 
do 
  case $yn in
      Yes ) echo "Write comment, then press ENTER:";
            read comment
            echo "$comment" >> comment.txt
            break;;
      No  ) echo "No comment." >> comment.txt
            break;;
  esac
done
# gives you the possibility to write a comment in a txt file

./OurCABA.EXT.x < ./$PARFILE

if [ $? -ne 0 ]
then
    echo "${RED}Some problem occurred while running OurCABA.${NORMAL}"
else
    echo "${GREEN}CABA code complete;${NORMAL}"
fi

sleep 1

# rm AdiabaticEnergies.dat 
# rm DiabaticEnergies.dat
# rm Pmat.dat 
# rm VQmat.dat 
# rm QuickEffV.dat 
# rm QuickPMat.dat  
# rm ParAndThresh.dat
# removes AtomIon1D output from working directory

check=0
i=1
while [ $check == 0 ]
do
    CHECKDIR=../../../../data_local/$AIDIR/${PARITY}${OPENCHANNELS}_${USED_BOUNDSTATES}bs_set${i}
    if [ -d "$CHECKDIR" ]
    then
        i=$((i+1))
    else
        check=1
    fi
done

DATADIR=${PARITY}${OPENCHANNELS}_${USED_BOUNDSTATES}bs_set${i}

mkdir ../../../../data_local/$AIDIR/$DATADIR
# creates new subdirectory to avoid overwriting data

cp ./$PARFILE ../../../../data_local/$AIDIR/$DATADIR/
# copies parameter file in data subdirectory

# mv QuickVQMat.dat ../../../../data_local/$AIDIR/
# cp VQmat_check.dat ../../../../data_local/$AIDIR/$DATADIR/
cp comment.txt ../../../../data_local/$AIDIR/$DATADIR/
cp Kmatrix*.dat ../../../../data_local/$AIDIR/$DATADIR/
cp Observables*.dat ../../../../data_local/$AIDIR/$DATADIR/
cp Smatrix*.dat ../../../../data_local/$AIDIR/$DATADIR/
cp TimeDelay*.dat ../../../../data_local/$AIDIR/$DATADIR/
cp grid.used.data ../../../../data_local/$AIDIR/$DATADIR/
# moves all output in data subdirectory

if [ $? -ne 0 ]
then
    echo "${RED}Some problem occurred while moving the files.${NORMAL}"
    exit 1
else 
    echo "Output moved subdirectory ${BLUE}$AIDIR${NORMAL} in ${BLUE}data_local${NORMAL};"
    cd ..
    rm -rf $WORKDIR
    echo "Working directory deleted;"
    echo "Number of open channels: ${BLUE}$OPENCHANNELS${NORMAL};" 
    echo "Number of used bound states: ${BLUE}$USED_BOUNDSTATES${NORMAL};"
    echo "Number of data set: ${BLUE}$i${NORMAL}."
fi

