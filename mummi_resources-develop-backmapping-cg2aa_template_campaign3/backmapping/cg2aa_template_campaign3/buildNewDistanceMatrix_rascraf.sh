#!/bin/bash

## This script creates Ras+CRaf protein distance restraint files

# - modules and $GMX variables should be exported prior to running this script

########################## Script control parameters

if [ $# != 10 ]; then
  echo "Usage: $0 [ORIGSTRUCT] [HELIX5STRUCTNEW] [HELIX5ENDORIG] [HELIX5ENDNOW] [CRDSTRUCTNEW] [CRDCARTRIDGESTART] [CRDCARTRIDGEEND] [CRDCARTRIDGEREPLACE] [CUTOFFDIST] [OFFSETLEVEL]"
  echo "       ORIGSTRUCT = structure file pre-feedback"
  echo "       HELIX5STRUCTNEW = structure file to be used for new Ras G domain helix 5 distance matrix information"
  echo "                         If you don't need to update this, then send the name of ORIGSTRUCT and set other options appropriately"
  echo "       HELIX5ENDORIG = the AA residue number (internal) that was the original last residue of helix 5 (should be 172 for C3 and 171 for C4)"
  echo "       HELIX5ENDNOW = the new end of Ras G domain helix 5"
  echo "       CRDSTRUCTNEW = structure file to be used for new CRD loop distance matrix information"
  echo "                         If you don't need to update this, then send the name of ORIGSTRUCT and set other options appropriately"
  echo "       CRDCARTRIDGESTART = the first residue of the CRD cartridge (should be 279)"
  echo "       CRDCARTRIDGEEND = the last residue of the CRD cartridge (should be 295)"
  echo "       CRDCARTRIDGEREPLACE = boolean, 1 to replace the CRD cartridge with feedback version and 0 to use original"
  echo "       CUTOFFDIST = distance in nanometers beyond which distance restraints will not be applied"
  echo "       OFFSETLEVEL = an integer such that distance restraints will only include residues i,j for j>i+OFFSETLEVEL"
  echo "                     to get everything within the cutoff distance, use OFFSETLEVEL=-1"
  exit
fi

ORIGSTRUCT=$1
HELIX5STRUCTNEW=$2
HELIX5ENDORIG=$3
HELIX5ENDNOW=$4
CRDSTRUCTNEW=$5
CRDCARTRIDGESTART=$6
CRDCARTRIDGEEND=$7
CRDCARTRIDGEREPLACE=$8
CUTOFFDIST=$9
OFFSETLEVEL=${10}

# Output file naming -- do not change this
pnam=rascraf

########################## Verification of inputs

# verify the input structures exist
if [ ! -e $ORIGSTRUCT ]; then
  echo "Error: could not find ORIGSTRUCT file $ORIGSTRUCT"
  exit
fi
if ((HELIX5ENDNOW>HELIX5ENDORIG)) && [ ! -e $HELIX5STRUCTNEW ]; then
  echo "Error: Helix 5 has grown, but the HELIX5STRUCTNEW file $HELIX5STRUCTNEW could not be found"
  exit
fi
if ((CRDCARTRIDGEREPLACE)) && [ ! -e $CRDSTRUCTNEW ]; then
  echo "Error: the CRDCARTRIDGEREPLACE boolean is on, but the CRDSTRUCTNEW file $CRDSTRUCTNEW could not be found"
  exit
fi

########################## Generate index file

# Note: ACE186 and LYS187 are never included in the distance restraints involving the N-terminal end of the Raf RBD because 
#       the ACE2(res187) at the N-terminus of this region in the CG simulations (C1*) is always very floppy
# Note: TRP321 and SER322 and NMA323 are never included in the distance restraints involving the C-terminal end of the Raf CRD because
#       the TRP321 and NMA322 at the C-terminus of this region in the CG simulations (C1*) is always very floppy
# Note: residues 265-268 are the RBD/CRD linker and are never included in the distance restraint matrix

rm -f tmpindex.ndx index.ndx
{
  if (( HELIX5ENDNOW != HELIX5ENDORIG )); then
    # the Ras G domain and Raf RBD -- for the new H5 length
    echo "r1-${HELIX5ENDNOW} | r188-264 & !aH*"
    echo "r1-${HELIX5ENDNOW} | r188-264 & aHN HA"
    echo "r1-${HELIX5ENDNOW} | r188-264 & rTHR ILE & aHB"
    echo "r CGW & !aH*"
    echo "\"r_1-${HELIX5ENDNOW}_r_188-264_&_!H*\" | \"r_1-${HELIX5ENDNOW}_r_188-264_&_HN_HA\" | \"r_1-${HELIX5ENDNOW}_r_188-264_&_THR_ILE_&_HB\" | r GTP MG | \"CGW_&_!H*\""
    # a backbone subselection for further enforcement 
    # (don't include a tertiary carbon and all its partners; e.g., exclude HA if have CA, C, N, CB)
    echo "r1-${HELIX5ENDNOW} | r188-264 & a HN N CA CB C O"
  fi
  if (( HELIX5ENDNOW >= HELIX5ENDORIG )); then
    # the Ras G domain and Raf RBD -- for the original H5 length
    echo "r1-${HELIX5ENDORIG} | r188-264 & !aH*"
    echo "r1-${HELIX5ENDORIG} | r188-264 & aHN HA"
    echo "r1-${HELIX5ENDORIG} | r188-264 & rTHR ILE & aHB"
    echo "r CGW & !aH*"
    echo "\"r_1-${HELIX5ENDORIG}_r_188-264_&_!H*\" | \"r_1-${HELIX5ENDORIG}_r_188-264_&_HN_HA\" | \"r_1-${HELIX5ENDORIG}_r_188-264_&_THR_ILE_&_HB\" | r GTP MG | \"CGW_&_!H*\""
    # a backbone subselection for further enforcement 
    # (don't include a tertiary carbon and all its partners; e.g., exclude HA if have CA, C, N, CB)
    echo "r1-${HELIX5ENDORIG} | r188-264 & a HN N CA CB C O"
  fi
  # the CRD, we will parse it later for the cassette
  echo "r269-320 & !aH*"
  echo "r269-320 & aHN HA"
  echo "r269-320 & rTHR ILE & aHB"
  echo "\"r_269-320_&_!H*\" | \"r_269-320_&_HN_HA\" | \"r_269-320_&_THR_ILE_&_HB\" | r ZN2"
  # a backbone subselection for further enforcement 
  # (don't include a tertiary carbon and all its partners; e.g., exclude HA if have CA, C, N, CB)
  echo "r269-320 & a HN N CA CB C O"

  echo q
} | $GMX make_ndx -f $ORIGSTRUCT -o tmpindex.ndx
# sed -i can be unexpectedly slow, so create a new file instead
sed "s/&/AND/g;s/\*/star/g;s/\!/not/g" tmpindex.ndx > index.ndx
rm -f tmpindex.ndx

########################## Create variables for the groups that are in the index file

## Note that the _Retain files are never used. They were initially created to support a new idea, but that was abandonded. Lefthere in case go back to use it later.

groupRasRafNow="r_1-${HELIX5ENDNOW}_r_188-264_AND_notHstar_r_1-${HELIX5ENDNOW}_r_188-264_AND_HN_HA_r_1-${HELIX5ENDNOW}_r_188-264_AND_THR_ILE_AND_HB_GTP_MG_CGW_AND_notHstar"
groupRasRafOrig="r_1-${HELIX5ENDORIG}_r_188-264_AND_notHstar_r_1-${HELIX5ENDORIG}_r_188-264_AND_HN_HA_r_1-${HELIX5ENDORIG}_r_188-264_AND_THR_ILE_AND_HB_GTP_MG_CGW_AND_notHstar"
groupRasRafNowRetain="r_1-${HELIX5ENDNOW}_r_188-264_AND_HN_N_CA_CB_C_O"
groupRasRafOrigRetain="r_1-${HELIX5ENDORIG}_r_188-264_AND_HN_N_CA_CB_C_O"
groupCRD="r_269-320_AND_notHstar_r_269-320_AND_HN_HA_r_269-320_AND_THR_ILE_AND_HB_ZN2"
groupCRDRetain="r_269-320_AND_HN_N_CA_CB_C_O"

########################## Make a list of Protein atoms vs residues

# - this will be used together with $OFFSETLEVEL to remove i,j interactions for i and j residues deemed too close to restrain
rm -f output.atom.res.protein
cat protein_${pnam}.itp|sed "s/;/; /;s/\[/[ /;s/\]/ ]/" | awk 'BEGIN{p=0}{if($1=="[" && $2 == "atoms")p=1; if(p && $1=="[" && $2 != "atoms")exit; if(p && $1!="[" && $1!=";")print $0}'|awk '{print $1,$3}' > output.atom.res.protein

########################## 
########################## There is a loop here over many sections.
########################## The point is to do the same thing for x and xRetain groups
########################## 

for buildRound in 1 2; do

  ########################## Define the inputs for the multiple rounds

  case $buildRound in
    1)
      # note that we built the retain versions first so that there will not be any naming overwrites
      usegroupRasRafNow=$groupRasRafNowRetain
      usegroupRasRafOrig=$groupRasRafOrigRetain
      usegroupCRD=$groupCRDRetain
      ;;
    2)
      usegroupRasRafNow=$groupRasRafNow
      usegroupRasRafOrig=$groupRasRafOrig
      usegroupCRD=$groupCRD
      ;;
    *)
      echo "Error: code mistake in case statement for buildRound variable. Exiting"
      exit
      ;;
  esac

  ########################## Make the distance restraint component for Ras + RBD + GTP + MG + CGW

  rm -f tmp.posre_${pnam}_part1.itp 
  rm -f posre_orig.itp posre_new.itp
  if (( HELIX5ENDNOW <= HELIX5ENDORIG )); then
    ### the part we will normally use
    # when the Ras G domain helix 5 ends at or before the point at which it originally ended, then we just truncate the distance matrix pairs at the new end
    echo $usegroupRasRafNow | $GMX genrestr -disre -f $ORIGSTRUCT -n index.ndx -o posre_orig.itp -cutoff $CUTOFFDIST -disre_dist 0 
    # make the temporary posre file for this section
    awk -v o=$OFFSETLEVEL 'NR==FNR{r[$1]=$2} NR!=FNR{if(NF==0||$1==";"||$1=="["){next}; a=r[$1]; b=r[$2]; if(a>=b - o && a<=b + o)next; print $1,$2,10,$6,$7,$8,1000 " ; residue_pair",a,b}' output.atom.res.protein posre_orig.itp > tmp.posre_${pnam}_part1.itp
  
  else
    # when the Ras G domain helix 5 gets longer, then we create a mashup
    # the original part is created normally and does not need special parsing since the HVR beyond helix 5 does not have distance restaints on it
    echo $usegroupRasRafOrig | $GMX genrestr -disre -f $ORIGSTRUCT -n index.ndx -o posre_orig.itp -cutoff $CUTOFFDIST -disre_dist 0
    # the new part is created and then parsed specially to remove everything we do not need
    echo $usegroupRasRafNow | $GMX genrestr -disre -f $HELIX5STRUCTNEW -n index.ndx -o posre_new.itp -cutoff $CUTOFFDIST -disre_dist 0
    {
      awk -v o=$OFFSETLEVEL 'NR==FNR{r[$1]=$2} NR!=FNR{if(NF==0||$1==";"||$1=="["){next}; a=r[$1]; b=r[$2]; if(a>=b - o && a<=b + o)next; print $1,$2,10,$6,$7,$8,1000 " ; residue_pair",a,b}' output.atom.res.protein posre_orig.itp
      awk -v o=$OFFSETLEVEL -v eo=$HELIX5ENDORIG -v en=$HELIX5ENDNOW 'NR==FNR{r[$1]=$2} NR!=FNR{if(NF==0||$1==";"||$1=="["){next}; a=r[$1]; b=r[$2]; if(a>=b - o && a<=b + o)next; if((a<=eo || a>en) && (b<=eo || b>en))next; print $1,$2,10,$6,$7,$8,1000 " ; residue_pair_feedbackH5",a,b}' output.atom.res.protein posre_new.itp
    } > tmp.posre_${pnam}_part1.itp
  fi
  rm -f posre_orig.itp posre_new.itp

  ########################## Make the distance restraint component for CRD + ZN2
  
  rm -f tmp.posre_${pnam}_part2.itp
  rm -f posre_orig.itp posre_new.itp
  if ((!CRDCARTRIDGEREPLACE)); then
    # when the CRD cartridge does not get replaced, we just make the entire distance matrix like the original
    echo $usegroupCRD | $GMX genrestr -disre -f $ORIGSTRUCT -n index.ndx -o posre_orig.itp -cutoff $CUTOFFDIST -disre_dist 0
    # make the temporary posre file for this section 
    awk -v o=$OFFSETLEVEL 'NR==FNR{r[$1]=$2} NR!=FNR{if(NF==0||$1==";"||$1=="["){next}; a=r[$1]; b=r[$2]; if(a>=b - o && a<=b + o)next; print $1,$2,10,$6,$7,$8,1000 " ; residue_pair",a,b}' output.atom.res.protein posre_orig.itp > tmp.posre_${pnam}_part2.itp
  
  else
    # when the CRD cartridge gets replaced, then we create a mashup
    echo $usegroupCRD | $GMX genrestr -disre -f $ORIGSTRUCT -n index.ndx -o posre_orig.itp -cutoff $CUTOFFDIST -disre_dist 0
    echo $usegroupCRD | $GMX genrestr -disre -f $CRDSTRUCTNEW -n index.ndx -o posre_new.itp -cutoff $CUTOFFDIST -disre_dist 0
    {
      # remove pairs involving one or more residues from the cartridge for the original part
      awk -v o=$OFFSETLEVEL -v cstart=$CRDCARTRIDGESTART -v cend=$CRDCARTRIDGEEND 'NR==FNR{r[$1]=$2} NR!=FNR{if(NF==0||$1==";"||$1=="["){next}; a=r[$1]; b=r[$2]; if(a>=b - o && a<=b + o)next; if((a>=cstart && a<=cend) || (b>=cstart && b<=cend))next; print $1,$2,10,$6,$7,$8,1000 " ; residue_pair",a,b}' output.atom.res.protein posre_orig.itp
      # only include pairs involving one or more residues from the cartridge for the original part
      awk -v o=$OFFSETLEVEL -v cstart=$CRDCARTRIDGESTART -v cend=$CRDCARTRIDGEEND 'NR==FNR{r[$1]=$2} NR!=FNR{if(NF==0||$1==";"||$1=="["){next}; a=r[$1]; b=r[$2]; if(a>=b - o && a<=b + o)next; if((a<cstart || a>cend) && (b<cstart || b>cend))next; print $1,$2,10,$6,$7,$8,1000 " ; residue_pair_feedbackCRD",a,b}' output.atom.res.protein posre_new.itp
    } > tmp.posre_${pnam}_part2.itp
  
  fi
  rm -f posre_orig.itp posre_new.itp

  ########################## Join the parts to make the posre_${pnam}.itp file we will actually use

  rm -f posre_${pnam}.itp
  {
    echo "[ bonds ]"
    cat tmp.posre_${pnam}_part1.itp 
    cat tmp.posre_${pnam}_part2.itp 
  } > posre_${pnam}.itp
  rm -f cat tmp.posre_${pnam}_part1.itp at tmp.posre_${pnam}_part2.itp

  ########################## These files should be huge; a failure may lead to basically empty file. If fewer than 10 lines, exit

  MINLINES=10
  if (( $(cat posre_${pnam}.itp | wc -l) < MINLINES )); then
    echo "Error: posre_${pnam}.itp had fewer than $MINLINES lines, this must be an error. Exiting and renaming the posre_${pnam}.itp file to _error_posre_${pnam}.itp_error_"
    echo "Error: the error occured in buildRound=$buildRound"
    mv posre_${pnam}.itp _error_posre_${pnam}.itp_error_
    exit
  fi

  ########################## move the retain versions to a new name

  # note that we built the retain versions first so that there will not be any naming overwrites
  if ((buildRound==1)); then
    mv posre_${pnam}.itp posre_${pnam}_retain.itp
  fi

done

########################## 
########################## The loop over the buildRound variable is over
##########################

########################## Make the original veraion (not the retain version) stronger for the first step 

rm -f posre_${pnam}_stronger.itp
cat posre_${pnam}.itp | sed "s/ 1000 ; residue_pair/ 10000 ; residue_pair/" > posre_${pnam}_stronger.itp

########################## clean up

rm -f output.atom.res.protein
rm -f posre_orig.itp posre_new.itp
rm -f index.ndx

