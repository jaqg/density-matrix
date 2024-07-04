#!/bin/bash

DAL="quadru"
#DAL="dalton-inp-alfredo"
MOL="H2O-1 H2O-2 H2O-3 H2O-4 H2O-5 H2O-6 H2O-7 H2O-8"
#MOL="H2O-1"

METHOD="MCSCF"

#WRK=/home/quinyonero/pruebas/MCSCF/H2O-scan
WRK=$(pwd)
SCR=/scr/quinyonero/pruebas/MCSCF/H2O
#PRG=/home/quinyonero/Release_2020/dalton
#BIN=build-jose-dynamic
PRG=/home/quinyonero/dalnep
BIN=build_q64
#PRG=/util/metd/Release_2020/dalton
#BIN=build_gcc92

export BASDIR=$PRG/basis/
export WRKMEM=4000000000

test -d $SCR || mkdir -p $SCR
cp $PRG/$BIN/dalton.x $SCR

for dal in $DAL ; do
    for mol in $MOL ; do
        echo "Running ${mol} for ${dal}"
        echo
#
        cd $WRK
        cp $dal.dal   $SCR/DALTON.INP
        cp $mol.mol   $SCR/MOLECULE.INP
#
        cd $SCR
        ./dalton.x
#
        cp CART.INP    $WRK/$mol.xyz
        cp DALTON.OUT  $WRK/${mol}.out
        mv DALTON.OUT  ${dal}-${mol}.out

        #cp density-mat-SIRI.dat $WRK/${mol}-density-mat-SIRI.dat

        # Density matrix for CC calculation:
	#cp AO-density-mat.dat $WRK/${dal}-${mol}-AO-dens-mat.dat

        #cp CAS-data $WRK/CAS-data
#       cp MCSCF-D1-matrix $WRK/${mol}-$METHOD-D1-matrix
#       cp MCSCF-D2-matrix $WRK/${mol}-$METHOD-D2-matrix
#       cp oneint_MObasis $WRK/${mol}-$METHOD-oneintMObasis
#       cp twoint_MObasis $WRK/${mol}-$METHOD-twointMObasis
	cp SIRIFC       $WRK/${mol}-$METHOD-SIRIFC
        cp molecule.fmt $WRK/${mol}-$METHOD-integrals
#
    done
done

exit 0

