#!/bin/bash
#BSUB -J _sinceCG_
#BSUB -nnodes 1
## Summit
## #BSUB -P LRN005
## #BSUB -W 02:00
## Lassen
#BSUB -P cancer
#BSUB -W 05:00

########################## Computer cluster specific parameters

#MYCLUSTER=KODIAK
#MYCLUSTER=SUMMIT
MYCLUSTER=LASSEN

# load in any necessary modules for setting up environment
if [ "$MYCLUSTER" == "KODIAK" ]; then
  module load gcc/8.3.0
  module load python/3.6-anaconda-5.0.1
  source activate mdaenv
fi
if [ "$MYCLUSTER" == "SUMMIT" ]; then
  module load gcc/7.4.0
  module load python/3.7.0
  source /autofs/nccs-svm1_proj/lrn005/spack3/share/spack/setup-env.sh
  spack load -r py-scipy
  spack load -r py-mdanalysis-mummi@mda0.19.2_with_ddcmd
  spack load -r gromacs /pk4r4em   ## single precision gromacs (gmx)
  spack load -r gromacs /3lqpmww   ## double precision gromacs (gmx_d)

  ### Info from Harsh on May 29, 2020:
  # pk4r4em               gromacs@2019.6.mummifix ~double_precision ~mpi
  # ktbyuy4               gromacs@2019.6.mummifix ~double_precision +mpi
  # 3lqpmww               gromacs@2019.6.mummifix +double_precision ~mpi
  # tl2nez4               gromacs@2019.6.mummifix +double_precision +mpi
fi
if [ "$MYCLUSTER" == "LASSEN" ]; then
  module load gcc/7.3.1
  module load python/3.7.2
  source /usr/gapps/kras/spack6/share/spack/setup-env.sh
  #spack compiler add /usr/tce/packages/gcc/gcc-7.3.1/
  spack load -r py-scipy
  spack load -r py-tqdm
  spack load -r py-mdanalysis-mummi@mda_1.0.1_ddcmd
  spack load -r gromacs /dvwa5ff  ## single precision gromacs (gmx)
  spack load -r gromacs /aldrilm   ## double precision gromacs (gmx_d)
fi

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
if [ "$MYCLUSTER" == "SUMMIT" ] || [ "$MYCLUSTER" == "LASSEN" ]; then
  ## These are the executables Chris installed on summit:
  #export GMXD=/gpfs/alpine/proj-shared/lrn005/cneale/exec/GROMACS/exec/gromacs-2019.6_limitEMstep0.05nm_accelerateGrompp_checkRESTRBONDS/serial_gnu_double/bin/gmx_d
  #export GMX=/gpfs/alpine/proj-shared/lrn005/cneale/exec/GROMACS/exec/gromacs-2019.6_limitEMstep0.05nm_accelerateGrompp_checkRESTRBONDS/serial_gnu/bin/gmx
  ## These are the executabels from spack -- going to have to assume that the spack modules are loaded such that the paths are set up correctly
  export GMXD=gmx_d
  export GMX=gmx
fi

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

# This is for Campaign 1*/3 inputs. A different script is necessary for Campaign 1 topology and structure files

./mod_topandpdb_campaign3.sh $CGTOP topol.top $CGPDB cg.pdb

########################## check if there was feedback and prepare the distance restraint input files

./master_make_new_distance_restraint_files.sh 

########################## run the CG-to-AA conversion script (sinceCG.sh)

for ((i=1;i<=NUMATTEMPTS;i++)); do
  TIME_START=$(date +%s)
  if [ "$MYCLUSTER" == "KODIAK" ]; then
    ./sinceCG.sh topol.top cg.pdb
  else
    jsrun -n1 -a1 -c32 -g0 -b none ./sinceCG.sh topol.top cg.pdb
  fi
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
 
