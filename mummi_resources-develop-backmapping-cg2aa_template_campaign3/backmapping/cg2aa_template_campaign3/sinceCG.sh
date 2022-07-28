#!/bin/bash

## This script constructs an output all-atom Charmm model based on an input Martini Coarse-grained model

# - see the ______TODO______ file for work that remains to be done

########################## Acknowledgements

# This script (sinceCG.sh) is based on Tsjerk Wassenaar's initram.sh VERSTAG=20150906-13-TAW downloaded from Tsjerk's git
# Modifications by Chris Neale in March - May 2020
# backward.py (called by this script) is from the same version of Tsjerk Wassenaar's code (with minor modifications)
# Mapping/ files are from Tsjerk, plus additions from Cesar Lopez and Chris Neale

########################## Get command-line inputs and error if they are not valid

# Usage: ./sinceCG.sh [CG top file] [CG pdb file]
# - test inputs can be found in EXAMPLE/
# - modules and $GMX, $GMXD variables should be exported prior to running this script
# - any necessary modifications to CG .top or .pdb files should be made prior to calling this script
#   examples of that processing can be found in the dorun.sh script
# - a variety of normal and special variant .itp files must be created prior to running this script
#   they are setup by createProteinTopology/make_all.sh
#   - if you want to save space after everything is set up, you can delete the local copy of createProteinTopology/ but:
#     (a) keep a copy of it somewhere, since it is necessary if you are going to modify any included topology files
#     (b) since charmm36-jun2015_DOE_2016_Dec20_martinimod_March2020.ff -> createProteinTopology/charmm36-jun2015_DOE_2016_Dec20_martinimod_March2020.ff
#         you would need to remove that symbolic link and replace it with a copy of the actual force field directory

if [ "$#" -ne 2 ] || [ ! -e $1 ] || [ ! -e $2 ]; then
  echo "Error: $0 requires exactly 2 command line arguments: [CG top file] [CG pdb file] -- and both input files must exist. Exiting"
  exit
fi
CGTOP=$1
CGPDB=$2

########################## Script control parameters

# - command-line values provided to gmx mdrun to control parallelization
#   *** BE CAREFUL ABOUT THE -rdd value. It can't be too large, or there will be a failure of domain decomposition.
#       it also can not be too small or the 1.2 nm setup distance restraints (actually, type 10 bonds) will not be communicated 
#       if the starting positions are terrible. 3.5 nm seems like a decent compromise for this size of system
#export NP="-nt 36 -dd 3 3 4 -rdd 3.5 -ntomp 1"
export NP="-nt 18 -dd 3 3 2 -rdd 3 -ntomp 1"

# - if ALLOW_CIS_ACYL_ESTER_ON_FINAL=1 then unlikely but not impossible stereochemical conditions in acyl ester bonds will not be marked as failure on the final step
#   Chris recommends ALLOW_CIS_ACYL_ESTER_ON_FINAL=0 because 5-us AA simulations from charmm-gui do not develop these issues
export ALLOW_CIS_ACYL_ESTER_ON_FINAL=0

# - if ALLOW_SAPI_INTERNAL_RING_CHANGE_ON_FINAL=1 then changes in inositol ring geometry will not be marked as failure on the final step
#   Chris recommends ALLOW_SAPI_INTERNAL_RING_CHANGE_ON_FINAL=1 beacuse 2-us AA simulations from charmm-gui do frequently develop these changes
export ALLOW_SAPI_INTERNAL_RING_CHANGE_ON_FINAL=1

# - set VERBOSEMDSTRING="-v" to get verbose output during all mdrun commands (good for debugging); set VERBOSEMDSTRING="" otherwise
#export VERBOSEMDSTRING="-v"
export VERBOSEMDSTRING=""

# stop gromacs from creating backup files and make it simply overwrite files
export GMX_MAXBACKUP=-1

########################## Print error message and quit for invalid script control parameter combinations (parameters defined above)

# nothing here

########################## Function definitions

function wait_file_write_completed {
  # in cases where we are concerned that a file write may not complete, call this function with the argument of filename

  SECONDSWAITSLEEP=1

  if (( $(echo $1 | grep "^/" | wc -l) == 1 )); then
    # path provided from root already
    fdir=$(basedir $1)
    fnam=$1 
  else
    # relative path provided, so prepend the current directory from root
    fdir=$(pwd)
    fnam=${fdir}/$1
  fi

  for ((;;)); do
    if (( $(lsof +d $fdir | grep $fnam | wc -l) == 0 )); then
      break
    fi
    echo "IO_NOTE: lsof waiting on $fnam"
    sleep $SECONDSWAITSLEEP
  done
}

###

function wait_and_poll {
  # use after putting a gromacs run in the background to catch errors and hangs
  cpid=$1         # the backgrounded job id
  nam=$2          # name to identify this part of the code
  outfile=$3      # name of the output file to be polled
  successfile=$4  # name of the file whose existence indicates success -- typically a .gro file for an mdrun command

  # we need to put some gmx calls in the background and poll them or else gromacs can sometimes hang. This is the polling frequency
  SECONDSWAITSLEEP=5

  for ((;;)); do
    if test -n "$(find . -maxdepth 1 -name 'step*pdb' -print -quit)" || test -n "$(find . -maxdepth 1 -name 'dd_dump_err*pdb' -print -quit)"; then
      echo "Error: found gromacs step files or dd_dump_err files during $nam ... killing $cpid"
      kill -9 $cpid
      # Normally better to hang here than to exit if the process didn't die because this way we can diagnose the problem more easily
      # but in the larger campaign run it seems better to simply exit, so comment out the wait statement below
      #wait
      exit
    fi
    if (( $(grep "inf," $outfile | wc -l) )); then
      echo "Error: found infinite force in $outfile ... killing $cpid"
      kill -9 $cpid
      # Normally better to hang here than to exit if the process didn't die because this way we can diagnose the problem more easily
      # but in the larger campaign run it seems better to simply exit, so comment out the wait statement below
      #wait
      exit
    fi
    if (( $(grep " nan " $outfile | wc -l) )); then
      echo "Error: found energy value of nan in $outfile ... killing $cpid"
      kill -9 $cpid
      # Normally better to hang here than to exit if the process didn't die because this way we can diagnose the problem more easily
      # but in the larger campaign run it seems better to simply exit, so comment out the wait statement below
      #wait
      exit
    fi
    if [ -e $successfile ]; then
      break
    fi
    # exit the loop if the process has finished (some uncaught error)
    if ! ps -p $cpid > /dev/null; then
      break
    fi
    sleep $SECONDSWAITSLEEP
  done

  # A wait here may be important to ensure that the process has finished writing the file
  wait

  # Ensure that the output file is no longer being accessed
  wait_file_write_completed $successfile
}

########################## Timing info

TIME_START=$(date +%s)

########################## Run backward.py to get the backmapped AA coordinates

rm -f projected.gro 0-backward.gro backmapped.top backmapped.ndx
B="./backward.py -f $CGPDB -raw projected.gro -o 0-backward.gro -kick 0.05 -sol -p $CGTOP -po backmapped.top -n backmapped.ndx -from martini -to charmm36 -solname SOL"
echo $B; $B

# I can not figure out how to get this to map from CG ZN bead to AA atom ZA with resname ZN2
# - therefore, just add the change here 
sed -i "s/ZN      ZN/ZN2     ZN/" 0-backward.gro

expectedfile=0-backward.gro
if [ ! -e "$expectedfile" ]; then
  echo "Error: failure to generate ${expectedfile}. Check for upstream error messages"
  exit
fi

########################## First EM -- there are no nonbonded interactions on lipid-lipid (but there are on protein-protein, different from Tsjerk's version)

expectedfile=0-backward.gro

# loop over initial EM with (a) no distance restraints, (b) regular distance restraints, (c) extra strong distance restraints
# - distance restraints are enabled via type 10 bonds
# Note regarding the number of steps in EMb
#   - the procedure works with 500 steps here
#   - however, when distance restraints pull things together across long distances, there are fewer stereochemical issues if EMb is longer (e.g., 2500 steps)

export OMP_NUM_THREADS=1

for BASE in 1-EMa 1-EMb 1-EMc; do
  # run grompp
  rm -f ${BASE}.tpr
  G="$GMX grompp -f ${BASE}.mdp -c $expectedfile -n backmapped.ndx -p backmapped.top -o ${BASE} -maxwarn 2"
  echo $G; $G 

  # run mdrun in background and poll for errors or completion
  # - output .g96 here because it has more precision (less chance of overlapping atoms)
  rm -f ${BASE}.g96 
  M="$GMXD mdrun -deffnm $BASE $NP -c ${BASE}.g96 $VERBOSEMDSTRING"
  echo $M 
  $M > output.mdrun.$BASE 2>&1 &
  cpid=$!
  wait_and_poll $cpid $BASE output.mdrun.$BASE ${BASE}.g96

  # variables we will use in the next round
  G96=${BASE}.g96

  expectedfile=$G96
  if [ ! -e "$expectedfile" ]; then
    echo "Error: failure to generate ${expectedfile}. Check for upstream error messages"
    exit
  fi
done

########################## Second EM -- includes all non-bonded interactions, but still using the topologies with modified potentials

LASTBASE=$BASE
BASE=2-EM

# run grompp
rm -f ${BASE}.tpr
G="$GMX grompp -f ${BASE}.mdp -c $G96 -n backmapped.ndx -p backmapped.top -o ${BASE} -maxwarn 2"
echo $G; $G

# run mdrun in background and poll for errors or completion
rm -f ${BASE}.gro
M="$GMXD mdrun -deffnm $BASE $NP $VERBOSEMDSTRING"
echo $M 
$M > output.mdrun.$BASE 2>&1 &
cpid=$!
wait_and_poll $cpid $BASE output.mdrun.$BASE ${BASE}.gro

expectedfile=${BASE}.gro
if [ ! -e "$expectedfile" ]; then
  echo "Error: failure to generate ${expectedfile}. Check for upstream error messages"
  exit
fi

########################## MD steps with increasing timestep

# Loop over the MD segments
# - .mdp files are already hard-coded into this directory, so changes to these loop values will require other changes
for DELTA_T in 0.0002 0.0005 0.001 0.002; do
  LASTBASE=$BASE
  BASE=3-MD-$DELTA_T

  # run grompp
  rm -f ${BASE}.tpr
  G="$GMX grompp -f ${BASE}.mdp -c ${LASTBASE}.gro -p backmapped.top -o ${BASE} -maxwarn 2"
  echo $G; $G
  
  # run mdrun in background and poll for errors or completion
  rm -f ${BASE}.gro
  M="$GMX mdrun -deffnm $BASE $NP $VERBOSEMDSTRING"
  echo $M 
  $M > output.mdrun.$BASE 2>&1 &
  cpid=$!
  wait_and_poll $cpid $BASE output.mdrun.$BASE ${BASE}.gro

  expectedfile=${BASE}.gro
  if [ ! -e "$expectedfile" ]; then
    echo "Error: failure to generate ${expectedfile}. Check for upstream error messages"
    exit
  fi
done  ## DELTA_T for-loop
 
# create a constant output file name from this point 
# Next step is mdanalysis, which does not understand about PBC, so must make a pbc -mol version
rm -f aa.gro
echo 0 | $GMX trjconv -f ${BASE}.gro -s ${BASE}.tpr -o aa.gro -pbc mol

########################## check and potentially fix stereochemistry 

# evaluate stereochemistry (check #1 of 3 possible)
python mdanalysis_chiral_cistrans.py aa.gro > output.mdanalysis
didfinish=$(grep "^Stereochemical Evaluation Complete" output.mdanalysis |wc -l)
if ((didfinish==0)); then
  echo "Error: some problem inhibited mdanalysis_chiral_cistrans.py from completing and writing out the final message: Stereochemical Evaluation Complete, exiting"
  exit
fi
if [ -e modified.gro ]; then
  mv modified.gro modified.1.gro
fi

# detect stereochemical issues
numissues=$(cat output.mdanalysis | grep -v "^Stereochemical Evaluation Complete" |wc -l)
numproteincisomega=$(cat output.mdanalysis | grep "cis omega" |wc -l)
# exit here if there is an internal ring issue... it is unlikely that this will be solved by anthing moving forward
shouldbail=$(cat output.mdanalysis | grep "unexpected internal ring"|wc -l)
if ((shouldbail)); then
  echo "Error: there is no code to recover from a bad internal ring dihedral, so exit"
  exit
fi

# if there was a stereochemical issue, run some EM on the structural modification that will hopefully fix the issue
if ((numissues!=0)); then
  # we have a modified file in modified.1.gro that we should energy minimize

  # run grompp
  rm -f chiralfix1.tpr
  if ((numproteincisomega>0)); then
    # there was at least one protein cis omega, so strengthen the dihedal for this EM
    cat backmapped.top | sed -E "s/protein_(.*)_CN_MARTINI_BACKMAPPING_ENFORCE.itp/protein_\1_CN_MARTINI_BACKMAPPING_ENFORCE_specialProteinOmega.itp/" > pco_backmapped.top
    G="$GMX grompp -f fixup.mdp -c modified.1.gro -p pco_backmapped.top -o chiralfix1.tpr -maxwarn 2"
  else
    # no protein cis omega, so use the regular toplogy
    G="$GMX grompp -f fixup.mdp -c modified.1.gro -p backmapped.top -o chiralfix1.tpr -maxwarn 2"
  fi
  echo $G; $G

  # run mdrun in background and poll for errors or completion
  rm -f chiralfix1.gro
  M="$GMXD mdrun -deffnm chiralfix1 $NP $VERBOSEMDSTRING"
  echo $M
  $M > output.mdrun.chiralfix1 2>&1 &
  cpid=$!
  wait_and_poll $cpid chiralfix1 output.mdrun.chiralfix1 chiralfix1.gro

  expectedfile=chiralfix1.gro
  if [ ! -e "$expectedfile" ]; then
    echo "Error: failure to generate ${expectedfile}. Check for upstream error messages"
    exit
  fi

  # evaluate stereochemistry to see if the modifications fixed it (check #2 of 3 possible)
  # Next step is mdanalysis, which does not understand about PBC, so must make a pbc -mol version
  rm -f chiralfix1_pbcmol.gro
  echo 0 | $GMX trjconv -f chiralfix1.gro -s chiralfix1.tpr -o chiralfix1_pbcmol.gro -pbc mol
  python mdanalysis_chiral_cistrans.py chiralfix1_pbcmol.gro > output.mdanalysis.2
  didfinish=$(grep "^Stereochemical Evaluation Complete" output.mdanalysis.2 |wc -l)
  if ((didfinish==0)); then
    echo "Error: some problem inhibited mdanalysis_chiral_cistrans.py from completing and writing out the final message: Stereochemical Evaluation Complete, exiting"
    exit
  fi
  if [ -e modified.gro ]; then
    mv modified.gro modified.2.gro
  fi
  numissues=$(cat output.mdanalysis.2 | grep -v "^Stereochemical Evaluation Complete" |wc -l)
  if ((numissues!=0)); then
    echo "Error: found chirality or double bond cis/trans issues on the second test, exiting"
    exit
  fi
  # there was an initial stereochemical error, but we fixed it so use the fix to move to EM with standard topologies
  nextstepgrofile=chiralfix1_pbcmol.gro
else 
  # there was no initial stereochemical eror, so we don't need to fix anything before moving to EM with standard topologies
  nextstepgrofile=aa.gro
fi

########################## copy structure file for convenience

# this so that we can check if final.gro exists for seeing how many tests failed -- we don't get here if there are any suspected problems

rm -f final.gro
ln -s $nextstepgrofile final.gro

########################## Create normal topology file

# Remove all of the special parameters we used to enforce desired stereochemistry
# - It seems possible that things will crash without an EM after removing the _CN_MARTINI_BACKMAPPING_ENFORCE parameters, 
#   but testing indicates things are stable and an extra EM is not necessary at this stage
cat backmapped.top | sed "s/_CN_MARTINI_BACKMAPPING_ENFORCE//g" > final.top

########################## velocity Langevin dynamics with no distance restraints and use the standard topology

# run grompp
rm -f LANGEVIN.tpr
G="$GMX grompp -f langevin.mdp -c final.gro -p final.top -o LANGEVIN.tpr -maxwarn 2"
echo $G; $G

# run mdrun in background and poll for errors or completion
rm -f LANGEVIN.gro
M="$GMX mdrun -deffnm LANGEVIN $NP $VERBOSEMDSTRING"
echo $M
$M > output.mdrun.LANGEVIN 2>&1 &
cpid=$!
wait_and_poll $cpid LANGEVIN output.mdrun.LANGEVIN LANGEVIN.gro

expectedfile=LANGEVIN.gro
if [ ! -e "$expectedfile" ]; then
  echo "Error: failure to generate ${expectedfile}. Check for upstream error messages"
  exit
fi

########################## check (but quit on error instead of trying to fix) stereochemistry (check #3 of 3 possible)

# - This is the final test before generating success.gro, so some stereochemical deviations may be allowed, depending on settings

# Next step is mdanalysis, which does not understand about PBC, so must make a pbc -mol version
rm -f LANGEVIN_pbcmol.gro
echo 0 | $GMX trjconv -f LANGEVIN.gro -s LANGEVIN.tpr -o LANGEVIN_pbcmol.gro -pbc mol
python mdanalysis_chiral_cistrans.py LANGEVIN_pbcmol.gro > output.mdanalysis.3
didfinish=$(grep "^Stereochemical Evaluation Complete" output.mdanalysis.3 |wc -l)
if ((didfinish==0)); then
  echo "Error: some problem inhibited mdanalysis_chiral_cistrans.py from completing and writing out the final message: Stereochemical Evaluation Complete, exiting"
  exit
fi
if [ -e modified.gro ]; then
  mv modified.gro modified.3.gro
fi

# Depending on settings, we may allow cis acy ester bonds and/or changes in inositol ring configurations -- we never allow changes in cholesterol rings
# - Chris suggests ALLOW_CIS_ACYL_ESTER_ON_FINAL=0 and ALLOW_SAPI_INTERNAL_RING_CHANGE_ON_FINAL=1
if (( ALLOW_CIS_ACYL_ESTER_ON_FINAL )); then
  if (( ALLOW_SAPI_INTERNAL_RING_CHANGE_ON_FINAL )); then
    numissues=$(cat output.mdanalysis.3 | grep -v "^Stereochemical Evaluation Complete" | grep -v "Error: cis acyl ester" | grep -v "Error: unexpected internal inositol ring" | wc -l)
  else
    numissues=$(cat output.mdanalysis.3 | grep -v "^Stereochemical Evaluation Complete" | grep -v "Error: cis acyl ester" | wc -l)
  fi
else
  if (( ALLOW_SAPI_INTERNAL_RING_CHANGE_ON_FINAL )); then
    numissues=$(cat output.mdanalysis.3 | grep -v "^Stereochemical Evaluation Complete" | grep -v "Error: unexpected internal inositol ring" | wc -l)
  else
    numissues=$(cat output.mdanalysis.3 | grep -v "^Stereochemical Evaluation Complete" | wc -l)
  fi
fi
if ((numissues!=0)); then
  echo "Error: found chirality or double bond cis/trans or internal ring dihedral issues on the third test, exiting"
  exit
fi

########################## copy structure file for convenience + ensure PBC has been solved to make molecules whole

# if success.gro exists, then we believe the build was successful; otherwise, repeat this script again
# - repeated runs of this script seem to have uncorrelated or weakly correlated chances of stereochemical error or crash

if [ -e LANGEVIN_pbcmol.gro ]; then
  rm -f success.gro
  ln -s LANGEVIN_pbcmol.gro success.gro
fi

########################## Print timing information 

TIME_END=$(date +%s)

echo "Total time required for build: $((TIME_END-TIME_START)) seconds"
echo "Time required for each mdrun segment:"
grep Time: *.log|awk '{print $1,$4/60}' | column -t
echo "NOTE: the stereochemical analysis/fixing within MDAnalysis also take 30 seconds each,"
echo "      and there are at least 2 (can be up to 3) such evaluations per build."
echo "      Plus, grompp takes some time, as do the few geometric transformations."

########################## Print information about deviation from restrained distances

# This relies on modification in check.cpp that make gromacs output information only about type 10 bond distances (normally not even printed).
# Must provide any .tpr file that has the type 10 bond restraints included 
#   - i.e., .mdp file that made that .tpr file must have had define=-DPOSRES or define=-DPOSRES_STRONGER in it

for f in 3-MD-0.001 success; do
  rm -f output.restraint_satisfaction.${f}
  $GMX check -s1 3-MD-0.001.tpr -f ${f}.gro 2>&1 | grep "^Distance between atoms" | sed "s/,//" | awk '{d=$NF-$(NF-3); if(d<0)d*=-1; print d,$0}' |sort -g -k 1 > output.restraint_satisfaction.${f}
  maxunsatisfy=$(tail -n 1 output.restraint_satisfaction.${f} | awk '{print $1}') 
  echo "MAXIMUM_RESTRAINT_UNSATISFACTION for ${f}.gro is $maxunsatisfy nanometers"
done

