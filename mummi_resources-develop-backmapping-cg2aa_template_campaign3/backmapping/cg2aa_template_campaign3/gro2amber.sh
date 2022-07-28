# setup gromacs
export GMX=`which gmx`

# setup amber
#. /usr/gapps/mummi/lassen/amber18/amber.sh
# not needed for summit (alreaddy sourced in summit . env)

# setup gromacs variables that assist amber 
export GMXLIB="$(pwd)"
topdir=$(echo $GMX | sed "s@bin/gmx@share/gromacs/top/@")

# define inputs
gmxtop=gro2amber.top
gmxgro=gro2amber.gro

# create the full.top file (let gromacs unroll all define statements prior to calling parmed)
rm -f tmp.tpr full.top mdout.mdp
touch empty.mdp
gmx grompp -f empty.mdp -p $gmxtop -c $gmxgro -o tmp.tpr -pp full.top
rm -f mdout.mdp

# call parmed to make the amber inputs
rm -f this.prmtop this.inpcrd
parmed << EOF
gromber full.top $gmxgro topdir $topdir
HMassRepartition
outparm amber.prmtop amber.inpcrd
checkValidity
EOF
