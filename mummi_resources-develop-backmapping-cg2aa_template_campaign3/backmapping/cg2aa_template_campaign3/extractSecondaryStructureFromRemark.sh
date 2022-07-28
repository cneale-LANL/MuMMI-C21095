#!/bin/bash

# this script takes a PDB file as input and returns the location of the end of Ras helix 5 and the end of the N-terminal CRD loop strand and the the start of the C-terminal CRD loop strand.
# - if it is just ras (no raf) then the CRD returns are -1
# - expected input is controlled by the first (and only) command line arg, which should be a string: ras, ras4A, rascraf, ras4Acraf

# newInputFile is a new AA structure file from which some regions will be taken for new AA structure definitions
# this must be a PDB file with one remark matching this format: 
#   - first string is "REMARK", second string is "300", and the third string is the secondary structure
#   - in the secondary structure listing, "2" and "H" are helix and "E" is sheet
# For Ras + Craf, the C1* secondary structure string would be: CEEEEEEEESTTSS1111HH2222SCCCCSCCCCSEEEEEEEEEETTEEEEEEEEEECCSTTCSS11112222CSEEEEEEETTC1111HHHHHHHHHH2222TCSCCCEEEEEECTTCSCCCSC1111HHH2222TCCEEECCTTTTCT1111HHHHHHHHHHHHH2222CCCCCCCCCCCCCCCCCEEEEEETTTEEEEEECCTTCC1111HH2222TTTCCCCCEEEEEECCSSSCSEEEECTTSCCCCCTTSCEEEEECTTCSEEECCEEEEECCSCCCCTTTCSCCSEEEEETTTCCEECCCCCCCSCSEEECCC
# Modifying the above secondary structure string only by adding some parenthetical indicators yields: C(<--Start of Ras at 2ACE1)EEEEEEEESTTSS1111HH2222SCCCCSCCCCSEEEEEEEEEETTEEEEEEEEEECCSTTCSS11112222CSEEEEEEETTC1111HHHHHHHHHH2222TCSCCCEEEEEECTTCSCCCSC1111HHH2222TCCEEECCTTTTCT1111HHHHHHHHHHHHH2222CCCCCCCCCCCCC(<--End of Ras at 185CYF)C(<--Start of Raf RBD at 187ACE2)CCCEEEEEETTTEEEEEECCTTCC1111HH2222TTTCCCCCEEEEEECCSSSCSEEEECTTSCCCCCTTSCEEEEEC(<--Start of Raf linker at 265LEU)TTC(<--End of Raf linker at 268VAL)SEEECCEEEEE(<--End of regular first sheet in CRD loops at 279THR)CCSCCCCTTTCSCCSE(<--Start of regular second sheet in CRD loops at 295ASN)EEEETTTCCEECCCCCCCSCSEEECCC(<--End of Raf CRD at 322NMA)
# now (below) on its own line with periods (or for missing residue 1 a hash) to mark out gaps so that counting will start at 1 and be continuous with character position 
#CEEEEEEEESTTSS1111HH2222SCCCCSCCCCSEEEEEEEEEETTEEEEEEEEEECCSTTCSS11112222CSEEEEEEETTC1111HHHHHHHHHH2222TCSCCCEEEEEECTTCSCCCSC1111HHH2222TCCEEECCTTTTCT1111HHHHHHHHHHHHH2222CCCCCCCCCCCCC.CCCCEEEEEETTTEEEEEECCTTCC1111HH2222TTTCCCCCEEEEEECCSSSCSEEEECTTSCCCCCTTSCEEEEECTTCSEEECCEEEEECCSCCCCTTTCSCCSEEEEETTTCCEECCCCCCCSCSEEECCC

if [ $# != 2 ]; then
  echo "Usage: $0 [NEWINPUTFILE] [PROTEINTYPE]"
  echo "       NEWINPUTFILE = the PDB file used to control this part of the feedback"
  echo "       PROTEINTYPE = a string, either: ras, ras4A, rascraf, or ras4Acraf"
  exit
fi

newInputFile=$1
proteinType=$2

case $proteinType in
  "ras")
    EXPECTEDLENGTH=184    # number of characters expected from secondary structure string
    STARTGHELIX=160       # position to start scanning the G domain helix 5 (in characters along the string)
    STARTCRDSHEETONE=-1   # position to start scanning the first marked CRD sheet (in characters along the string)
    ENDCRDSHEETTWO=-1     # position to start (backward) scanning the second marked CRD sheet (in characters along the string)
    ;;
  "rascraf")
    EXPECTEDLENGTH=320    # number of characters expected from secondary structure string
    STARTGHELIX=160       # position to start scanning the G domain helix 5 (in characters along the string)
    STARTCRDSHEETONE=273  # position to start scanning the first marked CRD sheet (in characters along the string)
    ENDCRDSHEETTWO=297    # position to start (backward) scanning the second marked CRD sheet (in characters along the string)
    ;;
  "ras4A")
    EXPECTEDLENGTH=185    # number of characters expected from secondary structure string
    STARTGHELIX=160       # position to start scanning the G domain helix 5 (in characters along the string)
    STARTCRDSHEETONE=-1   # position to start scanning the first marked CRD sheet (in characters along the string)
    ENDCRDSHEETTWO=-1     # position to start (backward) scanning the second marked CRD sheet (in characters along the string)
    ;;
  "ras4Acraf")
    EXPECTEDLENGTH=321    # number of characters expected from secondary structure string
    STARTGHELIX=160       # position to start scanning the G domain helix 5 (in characters along the string)
    STARTCRDSHEETONE=274  # position to start scanning the first marked CRD sheet (in characters along the string)
    ENDCRDSHEETTWO=298    # position to start (backward) scanning the second marked CRD sheet (in characters along the string)
    ;;
  "*")
    echo "Error: unexpected value of proteinType ($proteinType) provided to extractSecondaryStructureFromRemark.sh ; exiting"
    exit
    ;; 
esac

## ensure structure file exists
if [ ! -e $newInputFile ]; then
  echo "Error: structure file passed as argument to this script does not exist"
  exit
fi

## get the secondary structure
secstructstring=$(awk '{if($1=="REMARK" && $2==300 && NF==3){print $3}}' $newInputFile)
numfound=$(echo $secstructstring | wc -l)
if ((numfound!=1)); then
  echo "Error: expected 1 properly formatted REMARK 300 lines but found $numfound"
  exit
fi

## sanity check on length of the secondary structure string
secstructlength=$(echo $secstructstring | awk '{print length($1)}')
if (( secstructlength != EXPECTEDLENGTH )); then
  echo "Error: expected secondary structure length of $EXPECTEDLENGTH but found $secstructlength"
  exit
fi

## find the end of G domain helix 5
ghelixendsat=$(echo $secstructstring | awk '{lastres=$STARTGHELIX-1; for(i='$STARTGHELIX';i<='$secstructlength';i++){s=substr($1,i,1); if(s=="H" || s=="2")lastres=i; else exit}} END{print lastres}')
let realghelixendsat=$ghelixendsat+1
echo "The last residue of the G domain H5 is character $ghelixendsat which is really residue $realghelixendsat"

## find the end of the first marked CRD sheet
crdsheetoneendsat=$(echo $secstructstring | awk '{lastres=$STARTCRDSHEETONE-1; for(i='$STARTCRDSHEETONE';i<='$secstructlength';i++){s=substr($1,i,1); if(s=="E")lastres=i; else exit}} END{print lastres}')
let realcrdsheetoneendsat=$crdsheetoneendsat+2
echo "The last residue of the first marked CRD sheet is character $crdsheetoneendsat which is really residue $realcrdsheetoneendsat"

## find the start of the second marked CRD sheet
crdsheettwostartsat=$(echo $secstructstring | awk '{firstres=$ENDCRDSHEETTWO+1; for(i='$ENDCRDSHEETTWO';i>0;i--){s=substr($1,i,1); if(s=="E")firstres=i; else exit}} END{print firstres}')
let realcrdsheettwostartsat=$crdsheettwostartsat+2
echo "The first residue of the second marked CRD sheet is character $crdsheettwostartsat which is really residue $realcrdsheettwostartsat"

