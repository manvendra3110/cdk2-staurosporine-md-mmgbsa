#!/usr/bin/env bash
set -euo pipefail
NAME=$1
SMI=$2
NS=$3
W="md/$NAME"
mkdir -p "$W"
cp mdp/*.mdp protein.pdb "$W/"
python make_lig.py "$SMI" "dock/${NAME}_pose.pdb" "$W/lig.sdf" | tee "$W/charge.txt"
CH=$(awk '{print $3}' "$W/charge.txt")
cd "$W"

obabel lig.sdf -O lig.mol2
acpype -i lig.mol2 -b LIG -c bcc -n "$CH" -a gaff2

gmx pdb2gmx -f protein.pdb -o prot.gro -p topol.top -water tip3p -ff amber99sb-ildn -ignh
sed -i '/forcefield.itp"/a #include "LIG.acpype/LIG_GMX.itp"' topol.top
echo "LIG   1" >> topol.top

NP=$(sed -n 2p prot.gro)
NL=$(sed -n 2p LIG.acpype/LIG_GMX.gro)
{
  echo complex
  echo $((NP+NL))
  sed -n "3,$((NP+2))p" prot.gro
  sed -n "3,$((NL+2))p" LIG.acpype/LIG_GMX.gro
  tail -n1 prot.gro
} > complex.gro

gmx editconf -f complex.gro -o box.gro -c -d 1.0 -bt dodecahedron
gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2
echo SOL | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 1
gmx mdrun -deffnm em
echo -e "1 | 13\nname 17 Protein_LIG\n14 | 16\nname 18 Water_and_ions\nq" | gmx make_ndx -f em.gro -o index.ndx

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
gmx mdrun -deffnm nvt
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr
gmx mdrun -deffnm npt
sed "s/NSTEPS/$((NS*500000))/" md.mdp > md_run.mdp
gmx grompp -f md_run.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr
gmx mdrun -deffnm md

echo -e "Protein\nSystem" | gmx trjconv -s md.tpr -f md.xtc -o md_c.xtc -pbc mol -center
echo -e "Backbone\nSystem" | gmx trjconv -s md.tpr -f md_c.xtc -o md_fit.xtc -fit rot+trans

echo -e "Backbone\nBackbone" | gmx rms -s md.tpr -f md_fit.xtc -o rmsd_bb.xvg -tu ns
echo -e "Backbone\nUNL" | gmx rms -s md.tpr -f md_fit.xtc -o rmsd_lig.xvg -tu ns
echo "C-alpha" | gmx rmsf -s md.tpr -f md_fit.xtc -o rmsf.xvg -res
echo "Protein" | gmx gyrate -s md.tpr -f md_fit.xtc -o rg.xvg

echo "DONE $NAME"
