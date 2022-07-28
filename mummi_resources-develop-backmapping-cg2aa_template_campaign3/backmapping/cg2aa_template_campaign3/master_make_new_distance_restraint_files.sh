#!/bin/bash

## this script controls creation of the new distance restraint file for the build
## there are no command line inputs. Instead, the output is controlled by the non-existance, or existance + content of the feedback AA files
##
## Expecations:
## - 0 or 1 files matching RAS_RAF_TERNARY_scfix-CYFpos-HVR_*.pdb 
##   - 0 for no G domain helix 5 feedback
##   - 1 for G domain helix 5 feedback based on that structure
## - 0 or 1 files matching RAS_RAF_TERNARY_scfix-CYFpos-CRD_*.pdb
##   - 0 for no Raf CRD feedback
##   - 1 for Raf CRD feedback based on that structure
## - Any of the above PDB files, if provided should:
##   - be in standard PDB format
##   - have the same order of atoms as the $ORIGSTRUCT (though we accept different names for residues or atoms as long as they correspond)
##   - have the secondary structure string listed at the top in a REMARK 300 field (see extractSecondaryStructureFromRemark.sh for more details)

### The variable below is the only place that needs to be updated in this entire script deck to start Ras helix 5 from a dsifferent length
GLOBALSTARTHELIX5ENDORIG=171

################# Make the Ras+CRaf distance matrix

PROTEINTYPE=rascraf
ORIGSTRUCT=origFranken_rascraf.pdb
HELIX5ENDORIG=$GLOBALSTARTHELIX5ENDORIG
CRDCARTRIDGESTART=279
CRDCARTRIDGEEND=295
CUTOFFDIST=1.2
OFFSETLEVEL=-1

# find out how many files there are for ras G domain helix 5
HELIX5STRUCTNEW=$(find . -maxdepth 1 -name "RAS_RAF_TERNARY_scfix-CYFpos-HVR_*.pdb")
numHELIX5STRUCTNEW=$(echo $HELIX5STRUCTNEW | wc -w)
case $numHELIX5STRUCTNEW in
  0)
    # there is no feedback for Ras G domain helix 5
    HELIX5ENDNOW=$HELIX5ENDORIG
    # send the original struct name, even though it will not be used
    HELIX5STRUCTNEW=$ORIGSTRUCT
    ;;
  1)
    # there is feedback for Ras G domain helix 5
    HELIX5ENDNOW=$(./extractSecondaryStructureFromRemark.sh $HELIX5STRUCTNEW $PROTEINTYPE |grep "The last residue of the G domain H5 is character" | awk '{print $NF}')
    ;;
  *)
    # mummi error
    echo "Error: there should only be 0 or 1 files that match the pattern RAS_RAF_TERNARY_scfix-CYFpos-HVR_\*.pdb"
    echo "instead, here are the matches: $HELIX5STRUCTNEW"
    echo "exiting"
    exit
    ;;
esac

# find out how many files there are for raf CRD
CRDSTRUCTNEW=$(find . -maxdepth 1 -name "RAS_RAF_TERNARY_scfix-CYFpos-CRD_*.pdb")
numCRDSTRUCTNEW=$(echo $CRDSTRUCTNEW | wc -w)
case $numCRDSTRUCTNEW in
  0)
    # there is no feedback for Raf CRD
    CRDCARTRIDGEREPLACE=0
    # send the original struct name, even though it will not be used
    CRDSTRUCTNEW=$ORIGSTRUCT
    ;;
  1)
    # there is feedback for Raf CRD
    CRDCARTRIDGEREPLACE=1
    ;;
  *)
    # mummi error
    echo "Error: there should only be 0 or 1 files that match the pattern RAS_RAF_TERNARY_scfix-CYFpos-CRD_\*.pdb"
    echo "instead, here are the matches: $CRDSTRUCTNEW"
    echo "exiting"
    exit
    ;;
esac

# run the build procedure
CMD="./buildNewDistanceMatrix_rascraf.sh $ORIGSTRUCT $HELIX5STRUCTNEW $HELIX5ENDORIG $HELIX5ENDNOW $CRDSTRUCTNEW $CRDCARTRIDGESTART $CRDCARTRIDGEEND $CRDCARTRIDGEREPLACE $CUTOFFDIST $OFFSETLEVEL"
echo "Trying: $CMD"
$CMD

################# Make the Ras-only distance matrix

PROTEINTYPE=ras
ORIGSTRUCT=origFranken_ras.pdb
HELIX5ENDORIG=$GLOBALSTARTHELIX5ENDORIG
CUTOFFDIST=1.2
OFFSETLEVEL=-1

# find out how many files there are for ras G domain helix 5
HELIX5STRUCTNEW=$(find . -maxdepth 1 -name "RAS_scfix-CYFpos-HVR_*.pdb")
numHELIX5STRUCTNEW=$(echo $HELIX5STRUCTNEW | wc -w)
case $numHELIX5STRUCTNEW in
  0)
    # there is no feedback for Ras G domain helix 5
    HELIX5ENDNOW=$HELIX5ENDORIG
    # send the original struct name, even though it will not be used
    HELIX5STRUCTNEW=$ORIGSTRUCT
    ;;
  1)
    # there is feedback for Ras G domain helix 5
    HELIX5ENDNOW=$(./extractSecondaryStructureFromRemark.sh $HELIX5STRUCTNEW $PROTEINTYPE |grep "The last residue of the G domain H5 is character" | awk '{print $NF}')
    ;;
  *)
    # mummi error
    echo "Error: there should only be 0 or 1 files that match the pattern RAS_scfix-CYFpos-HVR_\*.pdb"
    echo "instead, here are the matches: $HELIX5STRUCTNEW"
    echo "exiting"
    exit
    ;;
esac

# run the build procedure
CMD="./buildNewDistanceMatrix_ras.sh $ORIGSTRUCT $HELIX5STRUCTNEW $HELIX5ENDORIG $HELIX5ENDNOW $CUTOFFDIST $OFFSETLEVEL"
echo "Trying: $CMD"
$CMD

################# Make the Ras4A-only distance matrix

PROTEINTYPE=ras4A
ORIGSTRUCT=origFranken_ras4A.pdb
HELIX5ENDORIG=$GLOBALSTARTHELIX5ENDORIG
CUTOFFDIST=1.2
OFFSETLEVEL=-1

# find out how many files there are for ras G domain helix 5
HELIX5STRUCTNEW=$(find . -maxdepth 1 -name "RAS4A_scfix-CYFpos-HVR_*.pdb")
numHELIX5STRUCTNEW=$(echo $HELIX5STRUCTNEW | wc -w)
case $numHELIX5STRUCTNEW in
  0)
    # there is no feedback for Ras G domain helix 5
    HELIX5ENDNOW=$HELIX5ENDORIG
    # send the original struct name, even though it will not be used
    HELIX5STRUCTNEW=$ORIGSTRUCT
    ;;
  1)
    # there is feedback for Ras G domain helix 5
    HELIX5ENDNOW=$(./extractSecondaryStructureFromRemark.sh $HELIX5STRUCTNEW $PROTEINTYPE |grep "The last residue of the G domain H5 is character" | awk '{print $NF}')
    ;;
  *)
    # mummi error
    echo "Error: there should only be 0 or 1 files that match the pattern RAS4A_scfix-CYFpos-HVR_\*.pdb"
    echo "instead, here are the matches: $HELIX5STRUCTNEW"
    echo "exiting"
    exit
    ;;
esac

# run the build procedure
CMD="./buildNewDistanceMatrix_ras4A.sh $ORIGSTRUCT $HELIX5STRUCTNEW $HELIX5ENDORIG $HELIX5ENDNOW $CUTOFFDIST $OFFSETLEVEL"
echo "Trying: $CMD"
$CMD

################# Make the Ras4A+CRaf distance matrix

PROTEINTYPE=ras4Acraf
ORIGSTRUCT=origFranken_ras4Acraf.pdb
HELIX5ENDORIG=$GLOBALSTARTHELIX5ENDORIG
CRDCARTRIDGESTART=280
CRDCARTRIDGEEND=296
CUTOFFDIST=1.2
OFFSETLEVEL=-1

# find out how many files there are for ras G domain helix 5
HELIX5STRUCTNEW=$(find . -maxdepth 1 -name "RAS4A_RAF_scfix-CYFpos-HVR_*.pdb")
numHELIX5STRUCTNEW=$(echo $HELIX5STRUCTNEW | wc -w)
case $numHELIX5STRUCTNEW in
  0)
    # there is no feedback for Ras G domain helix 5
    HELIX5ENDNOW=$HELIX5ENDORIG
    # send the original struct name, even though it will not be used
    HELIX5STRUCTNEW=$ORIGSTRUCT
    ;;
  1)
    # there is feedback for Ras G domain helix 5
    HELIX5ENDNOW=$(./extractSecondaryStructureFromRemark.sh $HELIX5STRUCTNEW $PROTEINTYPE |grep "The last residue of the G domain H5 is character" | awk '{print $NF}')
    ;;
  *)
    # mummi error
    echo "Error: there should only be 0 or 1 files that match the pattern RAS4A_RAF_scfix-CYFpos-HVR_\*.pdb"
    echo "instead, here are the matches: $HELIX5STRUCTNEW"
    echo "exiting"
    exit
    ;;
esac

# find out how many files there are for raf CRD
CRDSTRUCTNEW=$(find . -maxdepth 1 -name "RAS4A_RAF_scfix-CYFpos-CRD_*.pdb")
numCRDSTRUCTNEW=$(echo $CRDSTRUCTNEW | wc -w)
case $numCRDSTRUCTNEW in
  0)
    # there is no feedback for Raf CRD
    CRDCARTRIDGEREPLACE=0
    # send the original struct name, even though it will not be used
    CRDSTRUCTNEW=$ORIGSTRUCT
    ;;
  1)
    # there is feedback for Raf CRD
    CRDCARTRIDGEREPLACE=1
    ;;
  *)
    # mummi error
    echo "Error: there should only be 0 or 1 files that match the pattern RAS4A_RAF_scfix-CYFpos-CRD_\*.pdb"
    echo "instead, here are the matches: $CRDSTRUCTNEW"
    echo "exiting"
    exit
    ;;
esac

# run the build procedure
CMD="./buildNewDistanceMatrix_ras4Acraf.sh $ORIGSTRUCT $HELIX5STRUCTNEW $HELIX5ENDORIG $HELIX5ENDNOW $CRDSTRUCTNEW $CRDCARTRIDGESTART $CRDCARTRIDGEEND $CRDCARTRIDGEREPLACE $CUTOFFDIST $OFFSETLEVEL"
echo "Trying: $CMD"
$CMD



