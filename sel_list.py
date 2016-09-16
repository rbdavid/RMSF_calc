#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python

import sys
import MDAnalysis

pdb_file = sys.argv[1]
important = 'protein or nucleic or resname A5 A3 U5 atp adp PHX MG'

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
sugar = "name C5' H5' H5'' C4' H4' O4' C1' H1' C3' H3' C2' O2' HO2'" + " C5* H50 H51 C4* H40 O4* C1* H10 C3* H30 O3* H3' C2* H20 O2* H2'"		# DOES NOT INCLUDE THE O5' atom (which I will include in the phosphate atom selection string...
sugar_5= "name HO5' O5' or " + sugar
sugar_3= sugar + " O3' HO3' "
base = 'name N9 C8 H8 H80 N7 C5 C6 N6 H60 H61 H62 N1 C2 H2 N3 C4 O6 H1 H21 H22 H6 H5 N4 H41 H42 C2 O2 O4 H3'	# selection string that will select all appropriate atoms for any of the nucleic residues...

a_phos = 'name O5* O2A O1A PA O3A'
b_phos = 'name PB O1B O2B O3B'
g_phos = 'name PG O1G O2G O3G'
inorg_phos = 'name P O1 H1 O2 H2 O3 O4'

selection_list = []

### BEGINNING THE SELECTION LIST WITH PROTEIN SELECTIONS:
with open('testing_selection_list.txt','w') as f:
	for i in range(nRes):
		temp_resname = u_important.residues[i].resname
		if temp_resname in pro.resnames:
			selection_list.append(u_important.residues[i].select_atoms(sel1))
			f.write('%03d    %3s   %2d   %s\n' %(i,temp_resname,u_important.residues[i].n_atoms,sel1))
			selection_list.append(u_important.residues[i].select_atoms(sel2))
			f.write('%03d    %3s   %2d   %s\n' %(i,temp_resname,u_important.residues[i].n_atoms,sel2))
			selection_list.append(u_important.residues[i].select_atoms(sel3))
			f.write('%03d    %3s   %2d   %s\n' %(i,temp_resname,u_important.residues[i].n_atoms,sel3))
	
		elif temp_resname in nucleic.resnames:
			selection_list.append(u_important.residues[i].select_atoms(base))
			f.write('%03d    %3s   %2d   %s\n' %(i,temp_resname,u_important.residues[i].n_atoms,base))
			if temp_resname in ['A5','U5','C5','G5']:
				selection_list.append(u_important.residues[i].select_atoms(sugar_5))
				f.write('%03d    %3s   %2d   %s\n' %(i,temp_resname,u_important.residues[i].n_atoms,sugar_5))
				continue
			elif temp_resname in ['A3','U3','C3','G3']:
				selection_list.append(u_important.residues[i].select_atoms(sugar_3))
				f.write('%03d    %3s   %2d   %s\n' %(i,temp_resname,u_important.residues[i].n_atoms,sugar_3))
			else:
				selection_list.append(u_important.residues[i].select_atoms(sugar))
				f.write('%03d    %3s   %2d   %s\n' %(i,temp_resname,u_important.residues[i].n_atoms,sugar))
			selection_list.append(u_important.select_atoms("(resid %s and name P OP1 OP2 O5') or bynum %s" %(u_important.residues[i].resid,u_important.residues[i-1][-1].index+1)))
			f.write('%03d    %3s   %2d   Phosphate\n' %(i,temp_resname,u_important.residues[i].n_atoms))
			
		elif temp_resname in triphos.resnames:
			if temp_resname in ['atp','adp']:
				selection_list.append(u_important.residues[i].select_atoms(base))
				f.write('%03d    %3s   %2d   %s\n' %(i,temp_resname,u_important.residues[i].n_atoms,base))
				selection_list.append(u_important.residues[i].select_atoms(sugar))
				f.write('%03d    %3s   %2d   %s\n' %(i,temp_resname,u_important.residues[i].n_atoms,sugar))
			if temp_resname == 'atp':
				selection_list.append(u_important.residues[i].select_atoms(a_phos))
				f.write('%03d    %3s   %2d   %s\n' %(i,temp_resname,u_important.residues[i].n_atoms,a_phos))
				selection_list.append(u_important.residues[i].select_atoms(b_phos))
				f.write('%03d    %3s   %2d   %s\n' %(i,temp_resname,u_important.residues[i].n_atoms,b_phos))
				selection_list.append(u_important.residues[i].select_atoms(g_phos))
				f.write('%03d    %3s   %2d   %s\n' %(i,temp_resname,u_important.residues[i].n_atoms,g_phos))
			if temp_resname == 'adp':
				selection_list.append(u_important.residues[i].select_atoms(a_phos))
				f.write('%03d    %3s   %2d   %s\n' %(i,temp_resname,u_important.residues[i].n_atoms,a_phos))
				selection_list.append(u_important.residues[i].select_atoms(b_phos))
				f.write('%03d    %3s   %2d   %s\n' %(i,temp_resname,u_important.residues[i].n_atoms,b_phos))
			if temp_resname == 'PHX':
				selection_list.append(u_important.residues[i].select_atoms(inorg_phos))
				f.write('%03d    %3s   %2d   %s\n' %(i,temp_resname,u_important.residues[i].n_atoms,inorg_phos))
	
		elif temp_resname in other.resnames:
			selection_list.append(u_important.residues[i])
			f.write('%03d    %3s   %2d   MG\n' %(i,temp_resname,u_important.residues[i].n_atoms))

print len(selection_list)
