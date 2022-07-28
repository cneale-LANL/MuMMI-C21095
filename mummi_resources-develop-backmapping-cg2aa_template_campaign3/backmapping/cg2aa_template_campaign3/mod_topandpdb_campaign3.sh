#!/bin/bash

### This script takes as input a .top and .pdb file from a Campaign 1 simulation and provides the output .top and .pdb files 
### after processing them so that they can be used as input to sinceCG.sh.
### All libraries should have been loaded and $GMX defined prior to calling this script

INPUT_TOP=$1
OUTPUT_TOP=$2
INPUT_PDB=$3
OUTPUT_PDB=$4

########################## modify the topology file 

# Make a modified version of input topology $INPUT_TOP and output it as $OUTPUT_TOP
# - Main purpose is to align CG residue and bead names with expectations for input to sinceCG.sh
# - template.top is provided in this directory and used as an input that sets up all the .itp files
# - water and ions are removed from the .top file, but will still exist in the .pdb file and will be backmapped

rm -f $OUTPUT_TOP
{
  cat template.top
  cat $INPUT_TOP | sed "s/\[/[ /;s/;/; /" | awk 'BEGIN{p=0} {if(p && $1!=";")print $0;if($1=="[" && $2=="molecules")p=1}'|sed "s/POPX/POPC/;s/CHOL/CHL1/;s/PAP6/SAPI/;s/DPSM/SSM /"|sed "s/^W /;W /;s/^NA+ /;NA+ /;s/^CL- /;CL- /;s/^NA /;NA /;s/^CL /;CL /"
} > $OUTPUT_TOP

########################## modify the structure file

# Make a modified version of input structure $INPUT_PDB and output it as $OUTPUT_PDB
# - change residue names to match expectations of sinceCG.sh
# - change the names of the farnesyl beads to match Cesar's map file
# - add extra beads to get residue-specific implementations of ACE* and NMA* residues that Cesar has combined

# change atom and residue names as required (except for ACE* and NMA*)
rm -f tmp.modname.atoms.1.pdb
cat $INPUT_PDB | grep "^ATOM" | sed "s/POPX/POPC/;s/PAP6/SAPI/;s/DPSM/ SSM/;s/CHOL/CHL1/" | sed "s/F1  CYF/SC2 CYF/;s/F2  CYF/SC3 CYF/;s/F3  CYF/SC4 CYF/;s/F4  CYF/SC5 CYF/" | sed "s/P1  CYP/SC2 CYP/;s/P2  CYP/SC3 CYP/;s/P3  CYP/SC4 CYP/;s/P4  CYP/SC5 CYP/" | sed 's/./ /22' |sed "s/MG+ ION/MG   MG/;s/NA+ ION/NA   NA/;s/CL- ION/CL   CL/;s/HOH HOH/CGW CGW/" > tmp.modname.atoms.1.pdb

# change the residue names and insert beads related to NMA and NMA1
rm -f tmp.modname.atoms.2.pdb
cat tmp.modname.atoms.1.pdb | sed "s/NMA /NMA_/" | awk '
  BEGIN{
    gotnma=0
    gotnma1=0
  }
  {
    if(substr($0,14,2)=="BB"){
      # get the current BB coordinates
      x=substr($0,31,8)
      y=substr($0,39,8)
      z=substr($0,47,8)
      # if the "NMA " (now "NMA_") residue is found, save a modified version of its BB coordinates that push them away from previous residue by 1/2 that vector
      if(substr($0,18,4)=="NMA_"){
        gotnma=1
        newnmastring=sprintf("ATOM      0  BB  NMA     0    %8.3f%8.3f%8.3f  1.00  0.00", x+(x-px)/1.0, y+(y-py)/1.0, z+(z-pz)/1.0)
      }
      # if the "NMA1" residue is found, save a modified version of its BB coordinates that push them away from previous residue by 1/2 that vector
      if(substr($0,18,4)=="NMA1"){
        gotnma1=1
        newnma1string=sprintf("ATOM      0  BB  NMA     0    %8.3f%8.3f%8.3f  1.00  0.00", x+(x-px)/1.0, y+(y-py)/1.0, z+(z-pz)/1.0)
      }
      # save the previous BB coordinates
      px=x
      py=y
      pz=z
    }
    if(gotnma && substr($0,18,4)!="NMA_"){
      print newnmastring
      gotnma=0
    }
    if(gotnma1 && substr($0,18,4)!="NMA1"){
      print newnma1string
      gotnma1=0
    }
    print $0
  }
  END{
    # this is only inportant if these occur at the very end of ths ATOM input
    if(gotnma){
      print newnmastring
    }
    if(gotnma1){
      print newnma1string
    }
  }
' | sed "s/NMA_/SER /;s/NMA1/ASP /" > tmp.modname.atoms.2.pdb

## thought about using editconf -resnr 1 to assign residue numbers at this stage, but NAM res 0 followed by ACE res 0 will get renumbered properly (at least in gmx 5.1.5)
## therefore, that step is not necessary here

# change the residue names and insert beads related to ACE1, ACE2, and ACE3
# - use tac to process things in the reverse order, and will then invert them
rm -f tmp.modname.atoms.3.pdb
tac tmp.modname.atoms.2.pdb | awk '
  BEGIN{
    gotace1=0
    gotace2=0
    gotace3=0
  }
  {
    if(substr($0,14,2)=="BB"){
      # get the current BB coordinates
      x=substr($0,31,8)
      y=substr($0,39,8)
      z=substr($0,47,8)
      # if the "ACE1" residue is found, save a modified version of its BB coordinates that push them away from previous residue by 1/2 that vector
      if(substr($0,18,4)=="ACE1"){
        gotace1=1
        newace1string=sprintf("ATOM      0  BB  ACE     0    %8.3f%8.3f%8.3f  1.00  0.00", x+(x-px)/1.0, y+(y-py)/1.0, z+(z-pz)/1.0)
      }
      if(substr($0,18,4)=="ACE2"){
        gotace2=1
        newace2string=sprintf("ATOM      0  BB  ACE     0    %8.3f%8.3f%8.3f  1.00  0.00", x+(x-px)/1.0, y+(y-py)/1.0, z+(z-pz)/1.0)
      }
      if(substr($0,18,4)=="ACE3"){
        gotace3=1
        newace3string=sprintf("ATOM      0  BB  ACE     0    %8.3f%8.3f%8.3f  1.00  0.00", x+(x-px)/1.0, y+(y-py)/1.0, z+(z-pz)/1.0)
      }
      # save the previous BB coordinates
      px=x
      py=y
      pz=z
    }
    if(gotace1 && substr($0,18,4)!="ACE1"){
      print newace1string
      gotace1=0
    }
    if(gotace2 && substr($0,18,4)!="ACE2"){
      print newace2string
      gotace2=0
    }
    if(gotace3 && substr($0,18,4)!="ACE3"){
      print newace3string
      gotace3=0
    }
    print $0
  }
  END{
    # this is only inportant if these occur at the very end of ths ATOM input
    if(gotace1){
      print newace1string
    }
    if(gotace2){
      print newace2string
    }
    if(gotace3){
      print newace3string
    }
  }
' | sed "s/ACE1/THR /;s/ACE2/LYS /;s/ACE3/PRO /" > tmp.modname.atoms.3.pdb

# now print out the version that has atoms as desired
# - must use tac (not cat) to invert the inverted atom order
rm -f cg_modified.pdb $OUTPUT_PDB
{
  cat $INPUT_PDB | awk '{if($1=="ATOM")exit; print $0}'
  tac tmp.modname.atoms.3.pdb
  echo "TER"
  echo "ENDMDL"
} > cg_modified.pdb
$GMX editconf -f cg_modified.pdb -o $OUTPUT_PDB -resnr 1


