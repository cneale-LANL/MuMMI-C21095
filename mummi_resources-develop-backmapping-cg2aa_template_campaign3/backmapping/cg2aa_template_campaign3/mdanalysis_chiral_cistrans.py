import MDAnalysis as mda
import sys
import numpy as np

# fixing chiral errors works great
FIXCHIRAL=1
# fixing cis peptide bonds frequently fails
FIXCISTRANS=1
# no way at the moment to fix ring problems

# If a dihedral angle internal to a ring differs from the expectation value by more than DV_ERROR_CUTOFF_* values
# suggested values are 50 for cholesterol and 70 for pip2
DV_ERROR_CUTOFF_CHOL=50
DV_ERROR_CUTOFF_SAPI=70

# created the lipid ring dihedral listings (used near the end of this file) like this:
# cat ../template/createProteinTopology/getCHOL_dihedrals/out.param | awk '{v=$5+180; if(v>360)v-=360.0; print "grps.append([\""$1"\", \""$2"\", \""$3"\", \""$4"\", "v"])"}'


uu = mda.Universe(str(sys.argv[1]))
for ts in uu.trajectory: 
  # nested selection speeds up the algorithm by ~30% -- but be careful here if you are adding new lipid types
  u = uu.select_atoms('not resname SOL')
  sel_protein = u.select_atoms('protein or resname CYP CYF ACE NMA')
  sel_protein_cyp = sel_protein.select_atoms('resname CYP')
  sel_protein_cyf = sel_protein.select_atoms('resname CYF')
  sel_protein_nogly_noace_nonma = sel_protein.select_atoms('not resname GLY ACE NMA')
  sel_thr = sel_protein.select_atoms('resname THR')
  sel_ile = sel_protein.select_atoms('resname ILE')
  sel_lipids_withglycerol = u.select_atoms('resname POPS POPC POPE DIPE PAPC PAPS SAPI')
  sel_lipids_withps = sel_lipids_withglycerol.select_atoms('resname POPS PAPS')
  sel_lipids_ssm = u.select_atoms('resname SSM')
  sel_lipids_chol = u.select_atoms('resname CHL1')
  sel_lipids_sapi = sel_lipids_withglycerol.select_atoms('resname SAPI')
  sel_lipids_witholeoyl = sel_lipids_withglycerol.select_atoms('resname POPS POPC POPE')
  sel_lipids_dipe = sel_lipids_withglycerol.select_atoms('resname DIPE')
  sel_lipids_witharachidoyl = sel_lipids_withglycerol.select_atoms('resname PAPC PAPS SAPI')
  sel_gtp = u.select_atoms('resname GTP')

  didmod=0
  suggestbailout=0

  ##################################################################### PROTEIN
  ##################################################################### PROTEIN
  ##################################################################### PROTEIN

  ##################################################################### Protein chirality
  ##################################################################### Protein chirality
  ##################################################################### Protein chirality

  ## to get the desired chirality, put the hydrogen at the back, then provide the atoms in the desired order, always moving clockwise
  ## note that the order will be different for L and D enantiomers, but they should always be written clockwise here
  ## Chris Neale chose to use this order because this is how the map files for backward work in [chiral]: H, central atom, 3x listed clockwise with H at back

  ## Protein Calpha chirality
  for res in sel_protein_nogly_noace_nonma.residues:
    dih = res.atoms.select_atoms('name CA', 'name N', 'name CB', 'name C')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        print("Error: chirality Calpha issue for residue %s %d -- %f degrees" % (res.resname, res.resnum, dih.dihedral.value()))
        if FIXCHIRAL:
          heavystr='name CA'
          hydrogenstr='name HA'
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5 
          didmod+=1

  ## Threonine sidechain chirality
  for res in sel_thr.residues:
    dih = res.atoms.select_atoms('name CB', 'name CA', 'name CG2', 'name OG1')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        print("Error: chirality sidechain issue for residue %s %d -- %f degrees" % (res.resname, res.resnum, dih.dihedral.value()))
        if FIXCHIRAL:
          heavystr='name CB'
          hydrogenstr='name HB'
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5
          didmod+=1

  ## Isoleucine sidechain chirality
  for res in sel_ile.residues:
    dih = res.atoms.select_atoms('name CB', 'name CA', 'name CG2', 'name CG1')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        print("Error: chirality sidechain issue for residue %s %d -- %f degrees" % (res.resname, res.resnum, dih.dihedral.value()))
        if FIXCHIRAL:
          heavystr='name CB'
          hydrogenstr='name HB'
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5
          didmod+=1

  ## GTP chirality
  for res in sel_gtp.residues:
    dih = res.atoms.select_atoms('name C1\'', 'name N9', 'name C2\'', 'name O4\'')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        heavystr='name C1\''
        hydrogenstr='name H1\''
        print("Error: chirality GTP issue for residue %s %d (position of %s)-- %f degrees" % (res.resname, res.resnum, hydrogenstr, dih.dihedral.value()))
        if FIXCHIRAL:
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5
          didmod+=1
    dih = res.atoms.select_atoms('name C2\'', 'name C3\'', 'name O2\'', 'name C1\'')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        heavystr='name C2\''
        hydrogenstr='name H2\'\''
        print("Error: chirality GTP issue for residue %s %d (position of %s)-- %f degrees" % (res.resname, res.resnum, hydrogenstr, dih.dihedral.value()))
        if FIXCHIRAL:
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5
          didmod+=1
    dih = res.atoms.select_atoms('name C3\'', 'name C4\'', 'name O3\'', 'name C2\'')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        heavystr='name C3\''
        hydrogenstr='name H3\'\''
        print("Error: chirality GTP issue for residue %s %d (position of %s)-- %f degrees" % (res.resname, res.resnum, hydrogenstr, dih.dihedral.value()))
        if FIXCHIRAL:
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5
          didmod+=1
    dih = res.atoms.select_atoms('name C4\'', 'name O4\'', 'name C3\'', 'name C5\'')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        heavystr='name C4\''
        hydrogenstr='name H4\''
        print("Error: chirality GTP issue for residue %s %d (position of %s)-- %f degrees" % (res.resname, res.resnum, hydrogenstr, dih.dihedral.value()))
        if FIXCHIRAL:
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5
          didmod+=1


  ##################################################################### Protein cis/trans
  ##################################################################### Protein cis/trans
  ##################################################################### Protein cis/trans

  ## Protein backbone peptide bonds 
  #    - avoid "dih = res.omega_selection()" -- it is slow because it goes back to the universe each time
  #    - no explicit knowledge of separation between multiple consecutive protein chains
  #      Therefore, rely on ACE at start and CYF/NMA at end
  for res in sel_protein.residues:
    nextres=sel_protein.select_atoms("resnum %d"%(res.resnum+1))
    if len(nextres)==0:
      # this is the last amino acid found, so there is no i+1 residue
      dih=None
    else:
      if res.resname=="CYF" or res.resname=="NMA" or nextres.residues[0].resname=="ACE":
        # this dihedral would jump across protein chains, which is something that we want to avoid since it is nonsensical
        dih=None
      else:
        dih = (res.atoms.select_atoms('name CA', 'name C') + nextres.atoms.select_atoms('name N', 'name CA'))
        if len(dih) !=4:
          dih=None
    if dih is not None:
      v=dih.dihedral.value()
      if v>-90 and v<90:
        print("Error: cis omega for residue %s %d -- %f degrees" % (res.resname, res.resnum, dih.dihedral.value()))
        if FIXCISTRANS:
          nextres=sel_protein.select_atoms("resnum %d"%(res.resnum+1))
          castr='name CA'
          cstr='name C'
          ostr='name O'
          nstr='name N'
          hnstr='name HN'
          pos_ca = res.atoms.select_atoms(castr).positions
          pos_c = res.atoms.select_atoms(cstr).positions
          pos_o = res.atoms.select_atoms(ostr).positions
          pos_n = nextres.atoms.select_atoms(nstr).positions
          pos_hn = nextres.atoms.select_atoms(hnstr).positions
          vect_cac = pos_c - pos_ca
          vect_can = pos_n - pos_ca
          vect_cross = np.cross(vect_cac, vect_can)
          unit_vect_cross = vect_cross / np.linalg.norm(vect_cross)
          vect_co = pos_o - pos_c
          len_vect_co = np.linalg.norm(vect_co)
          res.atoms.select_atoms(ostr).positions = pos_c + unit_vect_cross * len_vect_co * 1.0
          nextres.atoms.select_atoms(hnstr).positions = pos_n - unit_vect_cross * len_vect_co * 1.0
          didmod+=1

  ## Ras farnesyl tail 2x trans dihedrals
  ## - I don't expect the fix to actually work, but let's see
  for res in sel_protein_cyf.residues:
    dih = res.atoms.select_atoms('name CJ', 'name CK', 'name CL', 'name CN')
    if dih is not None:
      v=dih.dihedral.value()
      if v>-90 and v<90:
        print("Error: cis Ras farnesyl CJ-CK-CL-CN for residue %s %d -- %f degrees" % (res.resname, res.resnum, dih.dihedral.value()))
        if FIXCISTRANS:
          mainstr='name CK'
          movestr='name HK'
          pos_main = res.atoms.select_atoms(mainstr).positions
          pos_move = res.atoms.select_atoms(movestr).positions
          vect_main_to_move = pos_move - pos_main
          res.atoms.select_atoms(movestr).positions -= vect_main_to_move*1.5
          didmod+=1
    dih = res.atoms.select_atoms('name CO', 'name CP', 'name CQ', 'name CS')
    if dih is not None:
      v=dih.dihedral.value()
      if v>-90 and v<90:
        print("Error: cis Ras farnesyl CO-CP-CQ-CS for residue %s %d -- %f degrees" % (res.resname, res.resnum, dih.dihedral.value()))
        if FIXCISTRANS:
          mainstr='name CP'
          movestr='name HP'
          pos_main = res.atoms.select_atoms(mainstr).positions
          pos_move = res.atoms.select_atoms(movestr).positions
          vect_main_to_move = pos_move - pos_main
          res.atoms.select_atoms(movestr).positions -= vect_main_to_move*1.5
          didmod+=1

  ## palmitoylated cys acyl thioester should be trans
  for res in sel_protein_cyp.residues:
    dih = res.atoms.select_atoms('name CB', 'name SG', 'name C1', 'name C2')
    if dih is not None:
      v=dih.dihedral.value()
      if v>-90 and v<90:
        print("Error: cis acyl thioester on palmitoylated cysteine for residue %s %d -- %f degrees" % (res.resname, res.resnum, dih.dihedral.value()))
        if FIXCISTRANS:
          mainstr='name C1'
          movestr='name O1'
          pos_main = res.atoms.select_atoms(mainstr).positions
          pos_move = res.atoms.select_atoms(movestr).positions
          vect_main_to_move = pos_move - pos_main
          res.atoms.select_atoms(movestr).positions -= vect_main_to_move*1.5
          didmod+=1


  ##################################################################### LIPIDS
  ##################################################################### LIPIDS
  ##################################################################### LIPIDS


  ##################################################################### Lipid chirality
  ##################################################################### Lipid chirality
  ##################################################################### Lipid chirality

  ## Lipid glycerol chirality
  for res in sel_lipids_withglycerol.residues:
    dih = res.atoms.select_atoms('name C2', 'name O21', 'name C1', 'name C3')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        heavystr='name C2'
        hydrogenstr='name HS'
        print("Error: chirality lipid glycerol issue for residue %s %d (position of %s)-- %f degrees" % (res.resname, res.resnum, hydrogenstr, dih.dihedral.value()))
        if FIXCHIRAL:
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5
          didmod+=1

  ## PS headgroup chirality
  for res in sel_lipids_withps.residues:
    dih = res.atoms.select_atoms('name C12', 'name N', 'name C11', 'name C13')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        heavystr='name C12'
        hydrogenstr='name H12A'
        print("Error: chirality PS headgroup issue for residue %s %d (position of %s)-- %f degrees" % (res.resname, res.resnum, hydrogenstr, dih.dihedral.value()))
        if FIXCHIRAL:
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5
          didmod+=1

  ## SSM backbone chirality
  for res in sel_lipids_ssm.residues:
    dih = res.atoms.select_atoms('name C2S', 'name C3S', 'name NF', 'name C1S')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        heavystr='name C2S'
        hydrogenstr='name H2S'
        print("Error: chirality SSM backbone issue for residue %s %d (position of %s)-- %f degrees" % (res.resname, res.resnum, hydrogenstr, dih.dihedral.value()))
        if FIXCHIRAL:
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5
          didmod+=1
    dih = res.atoms.select_atoms('name C3S', 'name C2S', 'name C4S', 'name O3')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        heavystr='name C3S'
        hydrogenstr='name H3S'
        print("Error: chirality SSM backbone issue for residue %s %d (position of %s)-- %f degrees" % (res.resname, res.resnum, hydrogenstr, dih.dihedral.value()))
        if FIXCHIRAL:
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5
          didmod+=1

  ## Cholesterol chirality
  for res in sel_lipids_chol.residues:
    dih = res.atoms.select_atoms('name C3', 'name O3', 'name C2', 'name C4')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        heavystr='name C3'
        hydrogenstr='name H3'
        print("Error: chirality cholesterol issue for residue %s %d (position of %s)-- %f degrees" % (res.resname, res.resnum, hydrogenstr, dih.dihedral.value()))
        if FIXCHIRAL:
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5
          didmod+=1
    dih = res.atoms.select_atoms('name C10', 'name C9', 'name C1', 'name C5')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        heavystr='name C10'
        hydrogenstr='name C19'
        print("Error: chirality cholesterol issue for residue %s %d (position of %s)-- %f degrees" % (res.resname, res.resnum, hydrogenstr, dih.dihedral.value()))
        if FIXCHIRAL:
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5
          ### >>>>>>>>>>>>>>>> Must move hydrogens attached to C19 as well
          res.atoms.select_atoms('name H19A').positions -= vect_heavy_to_hydrogen*2.0
          res.atoms.select_atoms('name H19B').positions -= vect_heavy_to_hydrogen*2.0
          res.atoms.select_atoms('name H19C').positions -= vect_heavy_to_hydrogen*2.0
          didmod+=1
    dih = res.atoms.select_atoms('name C9', 'name C8', 'name C10', 'name C11')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        heavystr='name C9'
        hydrogenstr='name H9'
        print("Error: chirality cholesterol issue for residue %s %d (position of %s)-- %f degrees" % (res.resname, res.resnum, hydrogenstr, dih.dihedral.value()))
        if FIXCHIRAL:
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5
          didmod+=1
    dih = res.atoms.select_atoms('name C8', 'name C7', 'name C14', 'name C9')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        heavystr='name C8'
        hydrogenstr='name H8'
        print("Error: chirality cholesterol issue for residue %s %d (position of %s)-- %f degrees" % (res.resname, res.resnum, hydrogenstr, dih.dihedral.value()))
        if FIXCHIRAL:
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5
          didmod+=1
    dih = res.atoms.select_atoms('name C14', 'name C8', 'name C13', 'name C15')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        heavystr='name C14'
        hydrogenstr='name H14'
        print("Error: chirality cholesterol issue for residue %s %d (position of %s)-- %f degrees" % (res.resname, res.resnum, hydrogenstr, dih.dihedral.value()))
        if FIXCHIRAL:
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5
          didmod+=1
    dih = res.atoms.select_atoms('name C13', 'name C17', 'name C12', 'name C14')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        heavystr='name C13'
        hydrogenstr='name C18'
        print("Error: chirality cholesterol issue for residue %s %d (position of %s)-- %f degrees" % (res.resname, res.resnum, hydrogenstr, dih.dihedral.value()))
        if FIXCHIRAL:
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5
          ### >>>>>>>>>>>>>>>> Must move hydrogens attached to C18 as well
          res.atoms.select_atoms('name H18A').positions -= vect_heavy_to_hydrogen*2.0
          res.atoms.select_atoms('name H18B').positions -= vect_heavy_to_hydrogen*2.0
          res.atoms.select_atoms('name H18C').positions -= vect_heavy_to_hydrogen*2.0
          didmod+=1
    dih = res.atoms.select_atoms('name C17', 'name C13', 'name C20', 'name C16')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        heavystr='name C17'
        hydrogenstr='name H17'
        print("Error: chirality cholesterol issue for residue %s %d (position of %s)-- %f degrees" % (res.resname, res.resnum, hydrogenstr, dih.dihedral.value()))
        if FIXCHIRAL:
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5
          didmod+=1
    dih = res.atoms.select_atoms('name C20', 'name C21', 'name C17', 'name C22')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        heavystr='name C20'
        hydrogenstr='name H20'
        print("Error: chirality cholesterol issue for residue %s %d (position of %s)-- %f degrees" % (res.resname, res.resnum, hydrogenstr, dih.dihedral.value()))
        if FIXCHIRAL:
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5
          didmod+=1

  ## SAPI headgroup chirality
  for res in sel_lipids_sapi.residues:
    dih = res.atoms.select_atoms('name C11', 'name O12', 'name C16', 'name C12')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        heavystr='name C11'
        hydrogenstr='name H1'
        print("Error: chirality SAPI headgroup issue for residue %s %d (position of %s)-- %f degrees" % (res.resname, res.resnum, hydrogenstr, dih.dihedral.value()))
        if FIXCHIRAL:
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5
          didmod+=1
    dih = res.atoms.select_atoms('name C12', 'name C11', 'name C13', 'name O2')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        heavystr='name C12'
        hydrogenstr='name H2'
        print("Error: chirality SAPI headgroup issue for residue %s %d (position of %s)-- %f degrees" % (res.resname, res.resnum, hydrogenstr, dih.dihedral.value()))
        if FIXCHIRAL:
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5
          didmod+=1
    dih = res.atoms.select_atoms('name C13', 'name C12', 'name C14', 'name O3')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        heavystr='name C13'
        hydrogenstr='name H3'
        print("Error: chirality SAPI headgroup issue for residue %s %d (position of %s)-- %f degrees" % (res.resname, res.resnum, hydrogenstr, dih.dihedral.value()))
        if FIXCHIRAL:
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5
          didmod+=1
    dih = res.atoms.select_atoms('name C14', 'name C13', 'name O4', 'name C15')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        heavystr='name C14'
        hydrogenstr='name H4'
        print("Error: chirality SAPI headgroup issue for residue %s %d (position of %s)-- %f degrees" % (res.resname, res.resnum, hydrogenstr, dih.dihedral.value()))
        if FIXCHIRAL:
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5
          didmod+=1
    dih = res.atoms.select_atoms('name C15', 'name C14', 'name C16', 'name O5')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        heavystr='name C15'
        hydrogenstr='name H5'
        print("Error: chirality SAPI headgroup issue for residue %s %d (position of %s)-- %f degrees" % (res.resname, res.resnum, hydrogenstr, dih.dihedral.value()))
        if FIXCHIRAL:
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5
          didmod+=1
    dih = res.atoms.select_atoms('name C16', 'name C15', 'name O6', 'name C11')
    if dih is not None:
      v=dih.dihedral.value()
      # could use something safer here, like "if v>-10" (since expect v ~-35), but we're going to try and invert the problematic ones
      if v>0:
        heavystr='name C16'
        hydrogenstr='name H6'
        print("Error: chirality SAPI headgroup issue for residue %s %d (position of %s)-- %f degrees" % (res.resname, res.resnum, hydrogenstr, dih.dihedral.value()))
        if FIXCHIRAL:
          pos_heavy = res.atoms.select_atoms(heavystr).positions
          pos_hydrogen = res.atoms.select_atoms(hydrogenstr).positions
          vect_heavy_to_hydrogen = pos_hydrogen - pos_heavy
          res.atoms.select_atoms(hydrogenstr).positions -= vect_heavy_to_hydrogen*1.5
          didmod+=1



  ##################################################################### Lipid cis/trans
  ##################################################################### Lipid cis/trans
  ##################################################################### Lipid cis/trans

  ## acyl ester should be trans (SN2)
  for res in sel_lipids_withglycerol.residues:
    dih = res.atoms.select_atoms('name C2', 'name O21', 'name C21', 'name C22')
    if dih is not None:
      v=dih.dihedral.value()
      if v>-90 and v<90:
        print("Error: cis acyl ester SN2 for residue %s %d -- %f degrees" % (res.resname, res.resnum, dih.dihedral.value()))
        if FIXCISTRANS:
          mainstr='name C21'
          movestr='name O22'
          pos_main = res.atoms.select_atoms(mainstr).positions
          pos_move = res.atoms.select_atoms(movestr).positions
          vect_main_to_move = pos_move - pos_main
          res.atoms.select_atoms(movestr).positions -= vect_main_to_move*1.5
          didmod+=1
  ## acyl ester should be trans (SN1)
  for res in sel_lipids_withglycerol.residues:
    dih = res.atoms.select_atoms('name C3', 'name O31', 'name C31', 'name C32')
    if dih is not None:
      v=dih.dihedral.value()
      if v>-90 and v<90:
        print("Error: cis acyl ester SN1 for residue %s %d -- %f degrees" % (res.resname, res.resnum, dih.dihedral.value()))
        if FIXCISTRANS:
          mainstr='name C31'
          movestr='name O32'
          pos_main = res.atoms.select_atoms(mainstr).positions
          pos_move = res.atoms.select_atoms(movestr).positions
          vect_main_to_move = pos_move - pos_main
          res.atoms.select_atoms(movestr).positions -= vect_main_to_move*1.5
          didmod+=1

  ## acyl chain oleoyl should be cis
  for res in sel_lipids_witholeoyl.residues:
    dih = res.atoms.select_atoms('name C28', 'name C29', 'name C210', 'name C211')
    if dih is not None:
      v=dih.dihedral.value()
      if v<-90 or v>90:
        print("Error: trans oleoyl for residue %s %d -- %f degrees" % (res.resname, res.resnum, dih.dihedral.value()))
        if FIXCISTRANS:
          mainstr='name C210'
          movestr='name H101'
          pos_main = res.atoms.select_atoms(mainstr).positions
          pos_move = res.atoms.select_atoms(movestr).positions
          vect_main_to_move = pos_move - pos_main
          res.atoms.select_atoms(movestr).positions -= vect_main_to_move*1.5
          didmod+=1

  ## acyl chain in lineoyl should have two cis
  for res in sel_lipids_dipe.residues:
    dih = res.atoms.select_atoms('name C28', 'name C29', 'name C210', 'name C211')
    if dih is not None:
      v=dih.dihedral.value()
      if v<-90 or v>90:
        print("Error: trans sn-2 proximal double bond in lineoyl for residue %s %d -- %f degrees" % (res.resname, res.resnum, dih.dihedral.value()))
        if FIXCISTRANS:
          mainstr='name C210'
          movestr='name H10R'
          pos_main = res.atoms.select_atoms(mainstr).positions
          pos_move = res.atoms.select_atoms(movestr).positions
          vect_main_to_move = pos_move - pos_main
          res.atoms.select_atoms(movestr).positions -= vect_main_to_move*1.5
          didmod+=1
    dih = res.atoms.select_atoms('name C38', 'name C39', 'name C310', 'name C311')
    if dih is not None:
      v=dih.dihedral.value()
      if v<-90 or v>90:
        print("Error: trans sn-1 proximal double bond in lineoyl for residue %s %d -- %f degrees" % (res.resname, res.resnum, dih.dihedral.value()))
        if FIXCISTRANS:
          mainstr='name C310'
          movestr='name H10X'
          pos_main = res.atoms.select_atoms(mainstr).positions
          pos_move = res.atoms.select_atoms(movestr).positions
          vect_main_to_move = pos_move - pos_main
          res.atoms.select_atoms(movestr).positions -= vect_main_to_move*1.5
          didmod+=1
    dih = res.atoms.select_atoms('name C211', 'name C212', 'name C213', 'name C214')
    if dih is not None:
      v=dih.dihedral.value()
      if v<-90 or v>90:
        print("Error: trans sn-2 distal double bond in lineoyl for residue %s %d -- %f degrees" % (res.resname, res.resnum, dih.dihedral.value()))
        if FIXCISTRANS:
          mainstr='name C213'
          movestr='name H13R'
          pos_main = res.atoms.select_atoms(mainstr).positions
          pos_move = res.atoms.select_atoms(movestr).positions
          vect_main_to_move = pos_move - pos_main
          res.atoms.select_atoms(movestr).positions -= vect_main_to_move*1.5
          didmod+=1
    dih = res.atoms.select_atoms('name C311', 'name C312', 'name C313', 'name C314')
    if dih is not None:
      v=dih.dihedral.value()
      if v<-90 or v>90:
        print("Error: trans sn-1 distal double bond in lineoyl for residue %s %d -- %f degrees" % (res.resname, res.resnum, dih.dihedral.value()))
        if FIXCISTRANS:
          mainstr='name C313'
          movestr='name H13X'
          pos_main = res.atoms.select_atoms(mainstr).positions
          pos_move = res.atoms.select_atoms(movestr).positions
          vect_main_to_move = pos_move - pos_main
          res.atoms.select_atoms(movestr).positions -= vect_main_to_move*1.5
          didmod+=1

  ## acyl chain in arachidoyl should have four cis
  for res in sel_lipids_witharachidoyl.residues:
    dih = res.atoms.select_atoms('name C24', 'name C25', 'name C26', 'name C27')
    if dih is not None:
      v=dih.dihedral.value()
      if v<-90 or v>90:
        print("Error: trans sn-2 0th from proximal double bond in arachidoyl for residue %s %d -- %f degrees" % (res.resname, res.resnum, dih.dihedral.value()))
        if FIXCISTRANS:
          mainstr='name C26'
          movestr='name H6R'
          pos_main = res.atoms.select_atoms(mainstr).positions
          pos_move = res.atoms.select_atoms(movestr).positions
          vect_main_to_move = pos_move - pos_main
          res.atoms.select_atoms(movestr).positions -= vect_main_to_move*1.5
          didmod+=1
    dih = res.atoms.select_atoms('name C27', 'name C28', 'name C29', 'name C210')
    if dih is not None:
      v=dih.dihedral.value()
      if v<-90 or v>90:
        print("Error: trans sn-2 1st from proximal double bond in arachidoyl for residue %s %d -- %f degrees" % (res.resname, res.resnum, dih.dihedral.value()))
        if FIXCISTRANS:
          mainstr='name C29'
          movestr='name H9R'
          pos_main = res.atoms.select_atoms(mainstr).positions
          pos_move = res.atoms.select_atoms(movestr).positions
          vect_main_to_move = pos_move - pos_main
          res.atoms.select_atoms(movestr).positions -= vect_main_to_move*1.5
          didmod+=1
    dih = res.atoms.select_atoms('name C210', 'name C211', 'name C212', 'name C213')
    if dih is not None:
      v=dih.dihedral.value()
      if v<-90 or v>90:
        print("Error: trans sn-2 2nd from proximal double bond in arachidoyl for residue %s %d -- %f degrees" % (res.resname, res.resnum, dih.dihedral.value()))
        if FIXCISTRANS:
          mainstr='name C212'
          movestr='name H12R'
          pos_main = res.atoms.select_atoms(mainstr).positions
          pos_move = res.atoms.select_atoms(movestr).positions
          vect_main_to_move = pos_move - pos_main
          res.atoms.select_atoms(movestr).positions -= vect_main_to_move*1.5
          didmod+=1
    dih = res.atoms.select_atoms('name C213', 'name C214', 'name C215', 'name C216')
    if dih is not None:
      v=dih.dihedral.value()
      if v<-90 or v>90:
        print("Error: trans sn-2 3rd from proximal double bond in arachidoyl for residue %s %d -- %f degrees" % (res.resname, res.resnum, dih.dihedral.value()))
        if FIXCISTRANS:
          mainstr='name C215'
          movestr='name H15R'
          pos_main = res.atoms.select_atoms(mainstr).positions
          pos_move = res.atoms.select_atoms(movestr).positions
          vect_main_to_move = pos_move - pos_main
          res.atoms.select_atoms(movestr).positions -= vect_main_to_move*1.5
          didmod+=1

  ## SSM C-linked trans
  for res in sel_lipids_ssm.residues:
    dih = res.atoms.select_atoms('name C3S', 'name C4S', 'name C5S', 'name C6S')
    if dih is not None:
      v=dih.dihedral.value()
      if v>-90 and v<90:
        print("Error: cis SSM C-linked double bond for residue %s %d -- %f degrees" % (res.resname, res.resnum, dih.dihedral.value()))
        if FIXCISTRANS:
          mainstr='name C5S'
          movestr='name H5S'
          pos_main = res.atoms.select_atoms(mainstr).positions
          pos_move = res.atoms.select_atoms(movestr).positions
          vect_main_to_move = pos_move - pos_main
          res.atoms.select_atoms(movestr).positions -= vect_main_to_move*1.5
          didmod+=1
  ## SSM N-linked trans -- I am not sure that I should enforce this, but can remove later
  for res in sel_lipids_ssm.residues:
    dih = res.atoms.select_atoms('name HNF', 'name NF', 'name C1F', 'name OF')
    if dih is not None:
      v=dih.dihedral.value()
      if v>-90 and v<90:
        print("Error: cis SSM N-linked resonance bond for residue %s %d -- %f degrees" % (res.resname, res.resnum, dih.dihedral.value()))
        if FIXCISTRANS:
          mainstr='name C1F'
          movestr='name OF'
          pos_main = res.atoms.select_atoms(mainstr).positions
          pos_move = res.atoms.select_atoms(movestr).positions
          vect_main_to_move = pos_move - pos_main
          res.atoms.select_atoms(movestr).positions -= vect_main_to_move*1.5
          didmod+=1



  ##################################################################### Lipid ring internal dihedrals
  ##################################################################### Lipid ring internal dihedrals
  ##################################################################### Lipid ring internal dihedrals


  ## Cholesterol ring internal dihedrals

  grps=[]

  # created the cholesterol ring dihedral listings like this:
  # cat ../template/createProteinTopology/getCHOL_dihedrals/out.param | awk '{v=$5+180; if(v>360)v-=360.0; print "grps.append([\""$1"\", \""$2"\", \""$3"\", \""$4"\", "v"])"}'

  grps.append(["C1", "C2", "C3", "C4", 56.42])
  grps.append(["C2", "C3", "C4", "C5", 304.263])
  grps.append(["C3", "C4", "C5", "C10", 53.219])
  grps.append(["C4", "C5", "C10", "C1", 311.581])
  grps.append(["C5", "C10", "C1", "C2", 49.629])
  grps.append(["C10", "C1", "C2", "C3", 304.743])
  grps.append(["C5", "C6", "C7", "C8", 13.1])
  grps.append(["C6", "C7", "C8", "C9", 316.537])
  grps.append(["C7", "C8", "C9", "C10", 57.807])
  grps.append(["C8", "C9", "C10", "C5", 320.178])
  grps.append(["C9", "C10", "C5", "C6", 9.029])
  grps.append(["C10", "C5", "C6", "C7", 4.462])
  grps.append(["C8", "C14", "C13", "C12", 298.483])
  grps.append(["C14", "C13", "C12", "C11", 56.349])
  grps.append(["C13", "C12", "C11", "C9", 307.82])
  grps.append(["C12", "C11", "C9", "C8", 47.723])
  grps.append(["C11", "C9", "C8", "C14", 312.147])
  grps.append(["C9", "C8", "C14", "C13", 57.515])
  grps.append(["C13", "C14", "C15", "C16", 328.471])
  grps.append(["C14", "C15", "C16", "C17", 6.789])
  grps.append(["C15", "C16", "C17", "C13", 19.582])
  grps.append(["C16", "C17", "C13", "C14", 322.353])
  grps.append(["C17", "C13", "C14", "C15", 43.172])
  grps.append(["C3", "C4", "C5", "C6", 235.941])
  grps.append(["C4", "C5", "C6", "C7", 181.547])
  grps.append(["C4", "C5", "C10", "C9", 191.813])
  grps.append(["C2", "C1", "C10", "C9", 171.39])
  grps.append(["C1", "C10", "C9", "C8", 200.088])
  grps.append(["C1", "C10", "C9", "C11", 68.754])
  grps.append(["C1", "C10", "C5", "C6", 128.797])
  grps.append(["C10", "C9", "C8", "C14", 179.814])
  grps.append(["C10", "C9", "C11", "C12", 178.825])
  grps.append(["C6", "C7", "C8", "C14", 195.233])
  grps.append(["C7", "C8", "C14", "C15", 307.222])
  grps.append(["C7", "C8", "C14", "C13", 179.921])
  grps.append(["C8", "C14", "C15", "C16", 196.893])
  grps.append(["C8", "C14", "C13", "C17", 176.012])
  grps.append(["C11", "C12", "C13", "C17", 167.315])
  grps.append(["C12", "C13", "C17", "C16", 207.437])
  grps.append(["C12", "C13", "C14", "C15", 165.643])
  
  for i in range(len(grps)):
    for res in sel_lipids_chol.residues:
      dih = res.atoms.select_atoms('name %s'%grps[i][0],'name %s'%grps[i][1],'name %s'%grps[i][2],'name %s'%grps[i][3])
      if dih is not None:
        # get the dihedral value
        v=dih.dihedral.value()
        if v<0.0:
          v+=360.0
        if v>=360.0:
          v-=360.0
        # get the delta from the restraint
        dv=v-grps[i][4]
        if dv<0.0:
          dv*=-1.0
        if dv>180.0:
          dv=360.0-dv
        if dv>DV_ERROR_CUTOFF_CHOL:
          print("Error: unexpected internal cholesterol ring (%s-%s-%s-%s) dihedral value for residue %s %d -- %.2f degrees (expected %.2f) -- |dDihe| = %.2f" % (grps[i][0],grps[i][1],grps[i][2],grps[i][3], res.resname, res.resnum, dih.dihedral.value(),grps[i][4],dv))



  ## Inositol ring internal dihedrals

  grps=[]

  # created the inositol ring dihedral listings like this:
  # cat ../template/createProteinTopology/getSAPI_dihedrals/out.param | awk '{v=$5+180; if(v>360)v-=360.0; print "grps.append([\""$1"\", \""$2"\", \""$3"\", \""$4"\", "v"])"}'

  grps.append(["C11", "C12", "C13", "C14", 59.296])
  grps.append(["C12", "C13", "C14", "C15", 301.376])
  grps.append(["C13", "C14", "C15", "C16", 58.428])
  grps.append(["C14", "C15", "C16", "C11", 301.141])
  grps.append(["C15", "C16", "C11", "C12", 59.316])
  grps.append(["C16", "C11", "C12", "C13", 300.453])

  for i in range(len(grps)):
    for res in sel_lipids_sapi.residues:
      dih = res.atoms.select_atoms('name %s'%grps[i][0],'name %s'%grps[i][1],'name %s'%grps[i][2],'name %s'%grps[i][3])
      if dih is not None:
        # get the dihedral value
        v=dih.dihedral.value()
        if v<0.0:
          v+=360.0
        if v>=360.0:
          v-=360.0
        # get the delta from the restraint
        dv=v-grps[i][4]
        if dv<0.0:
          dv*=-1.0
        if dv>180.0:
          dv=360.0-dv
        if dv>DV_ERROR_CUTOFF_SAPI:
          print("Error: unexpected internal inositol ring (%s-%s-%s-%s) dihedral value for residue %s %d -- %.2f degrees (expected %.2f) -- |dDihe| = %.2f" % (grps[i][0],grps[i][1],grps[i][2],grps[i][3], res.resname, res.resnum, dih.dihedral.value(),grps[i][4],dv))


  ##################################################################### Final considerations
  ##################################################################### Final considerations
  ##################################################################### Final considerations


  ## OUTPUT a new set of coordinates if made any chirality-based changes
  if didmod > 0:
    uu.atoms.write("modified.gro")

  ## print a line that we expect to find in sinceCG.sh -- purpose is to let that script fail if there is an error here
  print("Stereochemical Evaluation Complete")
