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
from distance_functions import *

zeros = np.zeros
sqrt = np.sqrt
flush = sys.stdout.flush

# ----------------------------------------
# VARIABLE DECLARATION

config_file = sys.argv[1]	# Local or Global positon of the config file that holds all the values for the parameters
necessary_parameters = ['avg_pdb','pdb_file','traj_loc','start','end','system','Wrapped','rmsf_filename','selection_file'] ###
all_parameters = ['avg_pdb','pdb_file','traj_loc','start','end','system','Wrapped','rmsf_filename','selection_file','alignment','important','substrate','protein_selection','write_summary','summary_filename']

sugar = "name C5' H5' H5'' C4' H4' O4' C1' H1' C3' H3' C2' O2' HO2'" + " C5* H50 H51 C4* H40 O4* C1* H10 C3* H30 O3* H3' C2* H20 O2* H2'"		# DOES NOT INCLUDE THE O5' atom (which I will include in the phosphate atom selection string...; the atoms with * are found in triphosphates;
sugar_5= "name HO5' O5' or " + sugar
sugar_3= sugar + " O3' HO3' "
base = 'name N9 C8 H8 H80 N7 C5 C6 N6 H60 H61 H62 N1 C2 H2 N3 C4 O6 H1 H21 H22 H6 H5 N4 H41 H42 C2 O2 O4 H3'	# selection string that will select all appropriate atoms for any of the nucleic residues...
a_phos = 'name O5* O2A O1A PA O3A'
b_phos = 'name PB O1B O2B O3B'
g_phos = 'name PG O1G O2G O3G'
inorg_phos = 'name P O1 H1 O2 H2 O3 O4'

# ----------------------------------------
# SUBROUTINES:

def ffprint(string):		# Useful function to use when on a computer that has a large buffer memory
	print '%s' %(string)
	flush()

def config_parser(config_file):	# Function to take config file and create/fill the parameter dictionary 
	for i in range(len(necessary_parameters)):
		parameters[necessary_parameters[i]] = ''

	# SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
	parameters['alignment'] = 'protein'
	parameters['important'] = 'protein'
	parameters['substrate'] = 'protein'
	parameters['protein_selection'] = 'not name H*'
	parameters['write_summary'] = False 
	parameters['summary_filename'] = 'rmsf.summary'

	# GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
	execfile(config_file,parameters)
	for key, value in parameters.iteritems():
		if value == '':
			print '%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.'
			sys.exit()
def summary():
	with open('%s.summary' %(parameters['summary_file']),'w') as f:
		f.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
		f.write('To recreate this analysis, run this line:\n')
		for i in range(len(sys.argv)):
			f.write('%s ' %(sys.argv[i]))
		f.write('\nParameters used:\n')
		for i in all_parameters:
			f.write('%s = %s \n' %(i, parameters[i]))
		f.write('\n\n')

# ----------------------------------------
# MAIN PROGRAM:
# CREATING PARAMETER DICTIONARY
parameters = {}
config_parser(config_file)

# ----------------------------------------
# LOAD IN REFERENCE STRUCTURE
ffprint('Beginning to prep the avg and u universes')
avg = MDAnalysis.Universe(parameters['avg_pdb'])
avg_align = avg.select_atoms(parameters['alignment'])
avg_important = avg.select_atoms(parameters['important'])
avg_important.translate(-avg_align.center_of_mass())	# translate the important selection by the COM of the alignment selection
pos0 = avg_align.positions		# gather the average position
nRes = avg_important.n_residues

# ----------------------------------------
# INITIALIZING THE ANALYSIS UNIVERSE; CREATING THE NECESSARY ATOM SELECTIONS FOR ALIGNMENT AND SUBSEQUENT SELECTION CREATION
u = MDAnalysis.Universe('%s' %(parameters['pdb_file']))
u_align = u.select_atoms(parameters['alignment'])
u_important = u.select_atoms(parameters['important'])
if not parameters['Wrapped']:
	u_substrate = u.select_atoms(parameters['substrate'])
	u_substrate_res = u_substrate.n_residues	#len(u_substrate.residues)

if nRes != u_important.n_residues:
	ffprint('Number of residues in the average structure: %d, Number of residues in the trajectory: %d. These need to match...' %(nRes,u_important.n_residues))
	sys.exit()

pro = u.select_atoms('protein')
nucleic = u.select_atoms('nucleic or resname A5 A3 U5 U3 G5 G3 C5 C3')
triphos = u.select_atoms('resname atp adp PHX')
other = u.select_atoms('resname MG')

# ----------------------------------------
# CREATING THE ATOM SELECTIONS FOR ANALYSIS
selection_list = []
nAtoms = []
avg_pos = []
with open('%s' %(parameters['selection_file']),'w') as f:
	for i in range(nRes):
		temp_resname = u_important.residues[i].resname
		if temp_resname in pro.resnames:
			temp_sel = u_important.residues[i].select_atoms(parameters['protein_selection'])
			selection_list.append(temp_sel)
			nAtoms.append(temp_sel.n_atoms)
			avg_pos.append(avg_important.residues[i].select_atoms(parameters['protein_selection']).positions)
			f.write('%03d    %3s   %2d   %s\n' %(i,temp_resname,temp_sel.n_atoms,parameters['protein_selection']))

		elif temp_resname in nucleic.resnames:
			temp_sel = u_important.residues[i].select_atoms(base)
			selection_list.append(temp_sel)
			nAtoms.append(temp_sel.n_atoms)
			avg_pos.append(avg_important.residues[i].select_atoms(base).positions)
			f.write('%03d    %3s   %2d   %s\n' %(i,temp_resname,temp_sel.n_atoms,base))
			if temp_resname in ['A5','U5','C5','G5']:
				temp_sel = u_important.residues[i].select_atoms(sugar_5)
				selection_list.append(temp_sel)
				nAtoms.append(temp_sel.n_atoms)
				avg_pos.append(avg_important.residues[i].select_atoms(sugar_5).positions)
				f.write('%03d    %3s   %2d   %s\n' %(i,temp_resname,temp_sel.n_atoms,sugar_5))
				continue

			elif temp_resname in ['A3','U3','C3','G3']:
				temp_sel = u_important.residues[i].select_atoms(sugar_3)
				selection_list.append(temp_sel)
				nAtoms.append(temp_sel.n_atoms)
				avg_pos.append(avg_important.residues[i].select_atoms(sugar_3).positions)
				f.write('%03d    %3s   %2d   %s\n' %(i,temp_resname,temp_sel.n_atoms,sugar_3))

			else:
				temp_sel = u_important.residues[i].select_atoms(sugar)
				selection_list.append(temp_sel)
				nAtoms.append(temp_sel.n_atoms)
				avg_pos.append(avg_important.residues[i].select_atoms(sugar).positions)
				f.write('%03d    %3s   %2d   %s\n' %(i,temp_resname,temp_sel.n_atoms,sugar))

		elif temp_resname in other.resnames:
			selection_list.append(u_important.residues[i])
			nAtoms.append(u_important.residues[i].n_atoms)
			avg_pos.append(avg_important.residues[i].positions)
			f.write('%03d    %3s   %2d   all\n' %(i,temp_resname,temp_sel.n_atoms))

ffprint('Done creating all the important atom selections, saving the number of atoms in the selection, and saving the average positions of the selections')
nSel = len(selection_list)
# ----------------------------------------
# MSD ANALYSIS OF THE ATOM SELECTIONS; SUMMING MSD VALUES OVER ALL TIMESTEPS;
dist2 = zeros(nSel)
ffprint('Beginning trajectory analysis')
nSteps = 0
start = parameters['start']
while start <= parameters['end']:
	ffprint('Loading/Analyzing trajectory %s' %(start))
	u.load_new('%s/production.%s/production.%s.dcd' %(parameters['traj_loc'],start,start))		### I STILL WOULD LIKE TO MAKE THIS LINE MORE GENERAL
	nSteps += len(u.trajectory)
	for ts in u.trajectory:
		u_important.translate(-u_align.center_of_mass())

		# CALCULATIONS that are unnecessary if the trajectory is wrapped.
		if not parameters['Wrapped']:		# Test to see if the 'Wrapped' key is equal to False
			dims = u.dimensions[:3]		
			dims2 = dims/2.0

			for i in range(u_substrate_res):
				COM = u_substrate.residues[i].center_of_mass()
				t = wrapping(COM,dims,dims2)
				u_substrate.residues[i].atoms.translate(t)

		R, d = rotation_matrix(u_align.positions,pos0)
		u_important.rotate(R)
	
		for i in range(nSel):
			dist2[i] += MSD(selection_list[i].positions,avg_pos[i],nAtoms[i])
		
		if ts.frame%1000 == 0:
			ffprint('Finished analyzing frame %d in trajectory %d.' %(ts.frame, start))

	start += 1

ffprint('Finished trajectory analysis.')
dist2 /= nSteps
dist2 = sqrt(dist2)

# ----------------------------------------
# WRITING RMSF RESULTS OUT TO FILE
with open('%s.dat' %(parameters['rmsf_filename']),'w') as f:
	np.savetxt(f,dist2)

if parameters['write_summary']:		# Test if 'write_summary' key is equal to True
	summary()

