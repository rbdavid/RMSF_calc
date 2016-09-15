#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python

import sys
import MDAnalysis

pdb_file = sys.argv[1]

u = MDAnalysis.Universe(pdb_file)
u_important = u.select_atoms(important)
nRes = u_important.n_residues

pro = u.select_atoms('protein')
nucleic = u.select_atoms('nucleic or resname A5 A3 U5 U3 G5 G3 C5 C3')
triphos = u.select_atoms('resname atp adp PHX')
other = u.select_atoms('resname MG')

sel1 = 'not name H*'
sel2 = 'backbone'
sel3 = 'name CA'

# SIMPLE SELECTION STRINGS TO USE FOR THE DESIRED SELECTIONS
sugar = " name C5' H5' H5'' C4' H4' O4' C1' H1' C3' H3' C2' O2' HO2' "	# DOES NOT INCLUDE THE O5' atom (which I will include in the phosphate atom selection string...
sugar_5= " name HO5' O5' or " + sugar
sugar_3= sugar + " O3' HO3' "
base = ' name N9 C8 H8 N7 C5 C6 N6 H61 H62 N1 C2 H2 N3 C4 O6 H1 H21 H22 H6 H5 N4 H41 H42 C2 O2 O4 H3'	# selection string that will select all appropriate atoms for any of the nucleic residues...

selection_list = []

### BEGINNING THE SELECTION LIST WITH PROTEIN SELECTIONS:
for i in range(nRes):
	temp_resname = u_important.residues[i].resname
	if temp_resname in pro.resnames:
		selection_list.append(u_important.residues[i].select_atoms(sel1))
		selection_list.append(u_important.residues[i].select_atoms(sel2))
		selection_list.append(u_important.residues[i].select_atoms(sel3))
	elif temp_resname in nucleic.resnames:
		selection_list.append(u_important.residues[i].select_atoms(base))
		if temp_resname in ['A5','U5','C5','G5']:
			selection_list.append(u_important.residues[i].select_atoms(sugar_5))
			continue
		elif temp_resname in ['A3','U3','C3','G3']:
			selection_list.append(u_important.residues[i].select_atoms(sugar_3))
		else:
			selection_list.append(u_important.residues[i].select_atoms(sugar))
		selection_list.append(u_important.select_atoms("(resid %s and anme P OP1 OP2 O5') or bynum %s" %(u_important.residues[i].resid,u_important.residues[i-1][-1].index+1)))
	elif temp_resname in triphos.resnames:
		### ... Develop these atom selections...

	elif temp_resname in other.resnames:
		selection_list.append(u_important.residues[i])

### SEE WHAT THIS OUTPUTS;...

