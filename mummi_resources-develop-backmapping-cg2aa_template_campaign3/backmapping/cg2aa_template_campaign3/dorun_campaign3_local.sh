#!/bin/bash

########################## Computer cluster specific parameters

# MuMMI host is set in the Mummi machine specific env file
MYCLUSTER=$MUMMI_HOST

#MYCLUSTER=KODIAK
#MYCLUSTER=summit
#MYCLUSTER=lassen

# define variables with full path to gromacs excutables
# - This package tested with gromacs version 2019.6
# - GMX is single precision (faster)
#   GMXD is double precision (can be necessary for energy minimizaiton)
# - Three source code modifications are necessary:
#   (A) modification to src/gromacs/mdlib/minimize.cpp near line 2595 to limit ustep to a maximum of 0.05 nm
#   (B) modification to src/gromacs/gmxpreprocess/convparm.c near line 558 that includes "|| i == F_RESTRBONDS" to make grompp faster
#   (C) modification to src/gromacs/tools/check.cpp near line 265 to print out F_RESTRBONDS (and only F_RESTRBONDS) distances
if [ "$MYCLUSTER" == "KODIAK" ]; then
  export GMXD=/net/scratch4/cneale/exe/KODIAK/GROMACS/exec/gromacs-2019.6_limitEMstep0.05nm_accelerateGrompp_checkRESTRBONDS/serial_double_gnu/bin/gmx_d
  export GMX=/net/scratch4/cneale/exe/KODIAK/GROMACS/exec/gromacs-2019.6_limitEMstep0.05nm_accelerateGrompp_checkRESTRBONDS/serial_gnu/bin/gmx
fi
if [ "$MYCLUSTER" == "summit" ]; then
  ## These are the executables Chris installed
  #export GMXD=/gpfs/alpine/proj-shared/lrn005/cneale/exec/GROMACS/exec/gromacs-2019.6_limitEMstep0.05nm_accelerateGrompp_checkRESTRBONDS/serial_gnu_double/bin/gmx_d
  #export GMX=/gpfs/alpine/proj-shared/lrn005/cneale/exec/GROMACS/exec/gromacs-2019.6_limitEMstep0.05nm_accelerateGrompp_checkRESTRBONDS/serial_gnu/bin/gmx
  ## These are the executabels from spack -- going to have to assume that the spack modules are loaded such that the paths are set up correctly
  export GMXD=gmx_d
  export GMX=gmx
fi

if [ "$MYCLUSTER" == "lassen" ]; then
  ## These are the executables Chris installed
  #export GMXD=/usr/gapps/kras/spack6/opt/spack/linux-rhel7-power9le/gcc-7.3.1/gromacs-2019.3-77ejhaxezkiewemzw3ig7yljyeq3lfrf/bin/gmx_d
  #export GMX=/usr/gapps/kras/spack6/opt/spack/linux-rhel7-power9le/gcc-7.3.1/gromacs-2019.6-5bbfaidl755223odzubaqqpilmenboey/bin/gmx
  export GMXD=gmx_d
  export GMX=gmx 
  ## These are the executabels from spack -- going to have to assume that the spack modules are loaded such that the paths are set up correctly
  #export GMXD=gmx_d
  #export GMX=gmx
fi

echo $MYCLUSTER
echo "gmx $GMX"
echo "gmx_d $GMXD"

########################## Definitions
    
# these are the names of the input files that were added to the template script deck for the cg-to-aa mapping
CGTOP=my.input.top
CGPDB=my.input.pdb
CGTPR=my.input.tpr

 ########################## Fix PBC and column that HOH starts in from Xiaohua's .pdb files

echo 0 | $GMX trjconv -f $CGPDB -o my.input.fix.pdb -s $CGTPR -pbc mol
CGPDB=my.input.fix.pdb

########################## Script control parameters 

# - attempt the backmapping until success, this is the maximum number of attempts
NUMATTEMPTS=3

########################## modify the topology and structure files to format them for expectations of sinceCG.sh

# The is for Campaign 1*/3 inputs. A different script is necessary for Campaign 1 topology and structure files

./mod_topandpdb_campaign3.sh $CGTOP topol.top $CGPDB cg.pdb

########################## check if there was feedback and prepare the distance restraint input files

./master_make_new_distance_restraint_files.sh

########################## run the CG-to-AA conversion script (sinceCG.sh)

for ((i=1;i<=NUMATTEMPTS;i++)); do
  TIME_START=$(date +%s)
  #jsrun -n1 -a1 -c32 -g0 -b none ./sinceCG.sh topol.top cg.pdb
  ./sinceCG.sh topol.top cg.pdb
  TIME_END=$(date +%s)
  if [ -e success.gro ]; then
    # successful AA system construction
    echo "SCRIPT_CONTROL: Round $i was a success."
    echo "SCRIPT_CONTROL: TIME TO SUCCESS = $((TIME_END-TIME_START)) seconds"
    break
  else
    # there was some error, try again
    echo "SCRIPT_CONTROL: Round $i was a failure. Directory contents listed in tmp.failure.file.listing.$i"
    echo "SCRIPT_CONTROL: TIME TO FAILURE = $((TIME_END-TIME_START)) seconds"
    ls -ltr > tmp.failure.file.listing.$i
    rm -f step*pdb dd_dump_err*pdb
  fi
done
 
