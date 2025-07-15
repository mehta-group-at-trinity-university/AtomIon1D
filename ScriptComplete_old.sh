#!/bin/bash


BOUNDSTATES=1
PARITY=P
OPENCHANNELS=1


# ––– DATA DIRECTORY NAME ––– 
AIDIR=${BOUNDSTATES}bs_R200

# ––– NAME OF ATOM ION PARAMETER FILE ––– 
AIPAR=AtomIon1D_${BOUNDSTATES}bs.par

# ––– NAME CABA PARAMETER FILE ––– 
CABAPAR=OurCABA_${PARITY}${OPENCHANNELS}.inp


# —————————————————————————————————————————————————————
# things below this line should not need to be modified

NORMAL=$(tput sgr0)
RED=$(tput setaf 1)
BLUE=$(tput setaf 4)
GREEN=$(tput setaf 2)
# colors for printed text


CABADIR=${BOUNDSTATES}bs 
# subdirectory in CABA

make
# compiles AtomIon1D.f90

if [ $? -ne 0 ]
then
    echo "${RED}AtomIon1D does not compile.${NORMAL}"
    exit 1
fi

if [ ! -d "../data_local/$AIDIR" ]
then
    mkdir ../data_local/$AIDIR
else
    printf "\nSubdirectory %s already exists.\nOverwrite? [Type number, then press ENTER]\n" "${BLUE}$AIDIR${NORMAL}"
    select yn in "Yes" "No"; 
    do 
      case $yn in
          Yes ) echo "${BLUE}Overwriting...${NORMAL}"; 
                sleep 2; 
                break;;
          No  )  echo "${RED}Change name of AtomIon subdirectory.${NORMAL}"; 
                exit 1;;
      esac
    done
fi
# creates subdirectory if it does not exist
# asks you to overwrite if already exists.

./AtomIon1D.x < ./input_AtomIon/$AIPAR
# runs AtomIon1D.x

if [ $? -ne 0 ]
then
    echo "${RED}Some problems occurred while running AtomIon1D.${NORMAL}"
    exit 1
fi

echo "${GREEN}Radial code complete;${NORMAL}"

check=0
i=1
while [ $check == 0 ]
do
    CHECKDIR=./CABA/${CABADIR}/workdir${i}
    if [ -d "$CHECKDIR" ]
    then
        i=$((i+1))
    else
        check=1
    fi
done

WORKDIR=workdir${i}
mkdir ./CABA/$CABADIR/$WORKDIR
echo "Working directory ${BLUE}$i${NORMAL} created;"

cp DiabaticEnergies.dat CABA/$CABADIR/$WORKDIR
cp AdiabaticEnergies.dat CABA/$CABADIR/$WORKDIR
cp VQmat.dat CABA/$CABADIR/$WORKDIR
cp Pmat.dat CABA/$CABADIR/$WORKDIR
cp QuickEffV.dat CABA/$CABADIR/$WORKDIR
cp QuickPMat.dat CABA/$CABADIR/$WORKDIR
cp ParAndThresh.dat CABA/$CABADIR/$WORKDIR
# copies files to working directory

cp ./input_AtomIon/$AIPAR ../data_local/$AIDIR/
mv DiabaticEnergies.dat ../data_local/$AIDIR/
mv AdiabaticEnergies.dat ../data_local/$AIDIR/
mv VQmat.dat ../data_local/$AIDIR/
mv Pmat.dat ../data_local/$AIDIR/
mv QuickEffV.dat ../data_local/$AIDIR/
mv QuickPMat.dat ../data_local/$AIDIR/
mv QTilmat.dat ../data_local/$AIDIR/
mv QuickQtil.dat ../data_local/$AIDIR/
mv ParAndThresh.dat ../data_local/$AIDIR/
# moves all output files subdirectory in case 
# you want to abort before the CABA is done

if [ $? -ne 0 ]
then
    echo "${RED}Problem while copying files.${NORMAL}"
    exit 1
fi

echo "Output files copied and moved to subdirectories;"

sleep 1

cd CABA
# goes to CABA directory
pwd

sleep 1

make
# compiles OurCABA.EXT.f

if [ $? -ne 0 ]
then
    echo "${RED}OurCABA does not compile.${NORMAL}"
    exit 1
fi

cp Legendre.dat ./$CABADIR/$WORKDIR
cp OurCABA.EXT.x ./$CABADIR/$WORKDIR
# copies executable and Legendre file to subdirectory

echo "CABA executable copied to working directory ${BLUE}$i${NORMAL} in ${BLUE}$CABADIR${NORMAL};"

sleep 1

cd $CABADIR/$WORKDIR
# goes to subdirectory
pwd

sleep 1

./OurCABA.EXT.x < ../input_files/$CABAPAR
# runs MyCABA.EXT.x

if [ $? -ne 0 ]
then
    echo "${RED}Some problem occurred while running OurCABA.${NORMAL}"
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
rm ParAndThresh.dat 
# removes all data files from AtomIon1D

check=0
i=1
while [ $check == 0 ]
do
    CHECKDIR=../../../../data_local/$AIDIR/${PARITY}${OPENCHANNELS}_set${i}
    if [ -d "$CHECKDIR" ]
    then
        i=$((i+1))
    else
        check=1
    fi
done

DATADIR=${PARITY}${OPENCHANNELS}_set${i}
mkdir ../../../../data_local/$AIDIR/$DATADIR
# creates new subdirectory to avoid overwriting data

cp ../input_files/$CABAPAR ../../../../data_local/$AIDIR/$DATADIR/
# copies CABA parameter file in data subdirectory

mv Kmatrix*.dat ../../../../data_local/$AIDIR/$DATADIR/
mv Observables*.dat ../../../../data_local/$AIDIR/$DATADIR/
mv Smatrix*.dat ../../../../data_local/$AIDIR/$DATADIR/
mv TimeDelay*.dat ../../../../data_local/$AIDIR/$DATADIR/
mv grid.used.data ../../../../data_local/$AIDIR/$DATADIR/
# moves all output in data subdirectory

echo "No comment." >> comment.txt
mv comment.txt ../../../../data_local/$AIDIR/$DATADIR/
# creates and moves comment file that is read in Mathematica


if [ $? -ne 0 ]
then
    echo "${RED}Some problem occurred while moving the files.${NORMAL}"
    exit 1
else 
    echo "Output moved subdirectory ${BLUE}$AIDIR${NORMAL} in ${BLUE}data_local${NORMAL};"
    cd ..
    rm -rf $WORKDIR
    echo "Working directory deleted;"
    echo "Number of open channels is ${BLUE}$OPENCHANNELS${NORMAL}, set number is ${BLUE}$i${NORMAL}." 
fi
