#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
##!/mnt/lustre_fs/users/mjmcc/apps/python2.7/bin/python
# ----------------------------------------
# USAGE:

# ----------------------------------------
# PREAMBLE:

import sys
import numpy as np
import MDAnalysis
from MDAnalysis.analysis.align import *
from sel_list import *
from distance_functions import *

# ----------------------------------------
# VARIABLE DECLARATION

avg_loc = sys.argv[1]
pdb_file = sys.argv[2]
traj_loc = sys.argv[3]
start = int(sys.argv[4])
end = int(sys.argv[5])
system = sys.argv[6]

alignment = 'protein and name CA and (resid 20:25 or resid 50:55 or resid 73:75 or resid 90:94 or resid 112:116 or resid 142:147 or resid 165:169 or resid 190:194 or resid 214:218 or resid 236:240 or resid 253:258 or resid 303:307)'
important = '(protein or nucleic or resname A5 or resname A3 or resname U5 or resname atp or resname adp or resname PHX or resname MG)'

zeros = np.zeros
sqrt = np.sqrt
flush = sys.stdout.flush

nSel = len(sel)

# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
	print '%s' %(string)
	flush()

# ----------------------------------------
# MAIN PROGRAM:
ffprint('Beginning to prep the avg and u universes')
# LOAD IN REFERENCE STRUCTURE
avg = MDAnalysis.Universe('%s%03d.%03d.avg_structure.pdb' %(avg_loc,start,end))
avg_align = avg.select_atoms(alignment)
avg_important = avg.select_atoms(important)

avg_important.translate(-avg_align.center_of_mass())
pos0 = avg_align.positions

nRes = len(avg_important.residues)

nAtoms = ['']*nRes
avg_pos = ['']*nRes 
avg_COM = zeros((nRes,3))
res_out = open('%s.residues_list.dat' %(system),'w')
for i in range(nRes):
	temp_list = []
	pos_list = []
	temp_res = avg_important.residues[i]
	res_out.write('%s   %d   %d\n' %(temp_res.resname, temp_res.resid, temp_res.resid+167))
	avg_COM[i] = temp_res.center_of_mass()
	for j in range(nSel):
		pos_list.append(temp_res.select_atoms('%s' %(sel[j][1])).positions)
		temp_list.append(len(temp_res.select_atoms('%s' %(sel[j][1]))))
	
	avg_pos[i] = pos_list
	nAtoms[i] = temp_list

res_out.close()

ffprint('Initialized and filled array with average residue coords. Also, created list of atoms to be analyzed for RMSF')

# LOAD IN PDB OF SYSTEM OF INTEREST
u = MDAnalysis.Universe(pdb_file)
u_align = u.select_atoms(alignment)
u_important = u.select_atoms(important)
u_substrate = u.select_atoms('nucleic or resname A5 or resname A3 or resname U5 or resname atp or resname adp or resname PHX or resname MG')

u_substrate_res = len(u_substrate.residues)

if nRes != len(u_important.residues):
	ffprint('Number of residues in the average structure: %d, Number of residues in the trajectory: %d. Obviously %d != %d. Fix it.' %(nRes,len(u_important.residues),nRes,len(u_important.residues)))
	sys.exit()

dist2 = zeros((nRes,nSel+1))
ffprint('Beginning trajectory analysis')
nSteps = 0
while start <= end:
	ffprint('Loading/Analyzing trajectory %s' %(start))
	u.load_new('%s/production.%s/production.%s.dcd' %(traj_loc,start,start))
	nSteps += len(u.trajectory)
	for ts in u.trajectory:
		dims = u.dimensions[:3]
		u_important.translate(-u_align.center_of_mass())

		for i in range(u_substrate_res):
			COM = zeros(3)
			COM = u_substrate.residues[i].center_of_mass()
			t = wrapping(COM,dims)
			u_substrate.residues[i].atoms.translate(t)

		R, d = rotation_matrix(u_align.positions,pos0)
		u_important.rotate(R)
		
		for i in range(nRes):
			temp_res = u_important.residues[i]
			for j in range(nSel):
				temp_pos = temp_res.select_atoms('%s' %(sel[j][1])).positions
#				print MSD(temp_pos,avg_pos[i][j],nAtoms[i][j]),nAtoms[i][j]
				dist2[i][j] += MSD(temp_pos,avg_pos[i][j],nAtoms[i][j])
			dist2[i][-1] += MSD(temp_res.center_of_mass(),avg_COM[i],1)
		
		if ts.frame%1000 == 0:
			ffprint('Finished analyzing frame %d in trajectory %d.' %(ts.frame, start))

	start += 1

ffprint('Finished trajectory analysis.')

dist2 /= nSteps

out1 = open('%s.COM_msd.dat' %(system),'w')
for i in range(nRes):
	out1.write('%f  \n' %(dist2[i][-1]))
out1.close()

dist2 = sqrt(dist2)

out2 = open('%s.residues_rmsf.dat' %(system),'w')
for i in range(nRes):
	for j in range(nSel):
		out2.write('%f   ' %(dist2[i][j]))
	out2.write('\n')
out2.close()

