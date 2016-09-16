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

zeros = np.zeros
sqrt = np.sqrt
flush = sys.stdout.flush

# ----------------------------------------
# VARIABLE DECLARATION

config_file = sys.argv[1]	# Local or Global positon of the config file that holds all the values for the parameters
necessary_parameters = ['avg_pdb','pdb_file','traj_loc','start','end','system','Wrapped','rmsf_filename','selection_file'] ###...
all_parameters = ['avg_pdb','pdb_file','traj_loc','start','end','system','Wrapped','rmsf_filename','selection_file','alignment','important','substrate','write_summary']
nSel = len(sel)

# ----------------------------------------
# SUBROUTINES:

def ffprint(string):		# Useful function to use when on a computer that has a large buffer memory
	print '%s' %(string)
	flush()

def config_parser(config_file):	# Function to take config file and create/fill the parameter dictionary 
	parameters = {}
	for i in range(len(necessary_parameters)):
		parameters[necessary_parameters[i]] = ''

	# SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
	parameters['alignment'] = 'protein'
	parameters['important'] = 'protein'
	parameters['substrate'] = 'protein'
	parameters['write_summary'] = False 

	# GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
	execfile(config_file,parameters)
	for key, value in parameters.iteritems:
		if value == '':
			print '%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.'
			sys.exit()
def summary():
		### ...

# ----------------------------------------
# MAIN PROGRAM:
# CREATING PARAMETER DICTIONARY
config_parser(config_file)

# LOAD IN REFERENCE STRUCTURE
ffprint('Beginning to prep the avg and u universes')
avg = MDAnalysis.Universe(parameters['avg_pdb'])
avg_align = avg.select_atoms(parameters['alignment'])
avg_important = avg.select_atoms(parameters['important'])
avg_important.translate(-avg_align.center_of_mass())	# translate the important selection by the COM of the alignment selection
pos0 = avg_align.positions		# gather the average position
nRes = avg_important.n_residues

nAtoms = ['']*nRes
avg_pos = ['']*nRes 
res_out = open('%s.residues_list.dat' %(parameters['system']),'w')
for i in range(nRes):
	temp_list = []
	pos_list = []
	temp_res = avg_important.residues[i]
	res_out.write('%s   %d   %d\n' %(temp_res.resname, temp_res.resid, temp_res.resid+167)) ### This offset value is very specific to Dengue NS3... need to make a general parameter to be read in...
	for j in range(nSel):
		pos_list.append(temp_res.select_atoms('%s' %(sel[j][1])).positions)
		temp_list.append(temp_res.select_atoms('%s' %(sel[j][1])).n_atoms)
	
	avg_pos[i] = pos_list
	nAtoms[i] = temp_list

res_out.close()

ffprint('Initialized and filled array with average residue coords. Also, created list of atoms to be analyzed for RMSF')	### NEED TO REWORD THIS OUTPUT...

# LOAD IN PDB OF SYSTEM OF INTEREST
u = MDAnalysis.Universe(parameters['pdb_file'])
u_align = u.select_atoms(parameters['alignment'])
u_important = u.select_atoms(parameters['important'])
if not parameters['Wrapped']:
	u_substrate = u.select_atoms(parameters['substrate'])
	u_substrate_res = u_substrate.n_residues	#len(u_substrate.residues)

if nRes != u_important.n_residues:
	ffprint('Number of residues in the average structure: %d, Number of residues in the trajectory: %d. These need to match...' %(nRes,u_important.n_residues))
	sys.exit()

dist2 = zeros((nRes,nSel))
ffprint('Beginning trajectory analysis')
nSteps = 0
start = parameters['start']
while start <= parameters['end']:
	ffprint('Loading/Analyzing trajectory %s' %(start))
	u.load_new('%s/production.%s/production.%s.dcd' %(parameters['traj_loc'],start,start))
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
		
		for i in range(nRes):
			temp_res = u_important.residues[i]
			for j in range(nSel):
				temp_pos = temp_res.select_atoms('%s' %(sel[j][1])).positions
				dist2[i][j] += MSD(temp_pos,avg_pos[i][j],nAtoms[i][j])
		
		if ts.frame%1000 == 0:
			ffprint('Finished analyzing frame %d in trajectory %d.' %(ts.frame, start))

	start += 1

ffprint('Finished trajectory analysis.')
dist2 /= nSteps
dist2 = sqrt(dist2)
# WRITING RMSF RESULTS OUT TO FILE
with open('%s.dat' %(parameters['rmsf_filename']),'w') as f:
	np.savetxt(f,dist2)

