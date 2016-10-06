#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:

# ----------------------------------------
# PREAMBLE:

from plotting_functions import *
import sys
import os 
import MDAnalysis

start = int(sys.argv[1])
end = int(sys.argv[2])
step = int(sys.argv[3])
pdb_file = sys.argv[4]
selection_file = sys.argv[5]
system = sys.argv[6]

change_dir = os.chdir

# ----------------------------------------
# MAIN PROGRAM

u = MDAnalysis.Universe(pdb_file)
pro = u.select_atoms('protein')
rna = u.select_atoms('nucleic or resname A5 U5 A3')
other = u.select_atoms('resname atp adp PHX MG')

selections = open(selection_file,'r')

pro_index = []
rna_index = []
other_index = []
nSelections = 0
for line in selections:
	line = line.split()
	temp_resname = line[1]
	nSelections += 1
	if temp_resname in pro.resnames:
		pro_index.append(int(line[0]))
	elif temp_resname in rna.resnames:
		rna_index.append(int(line[0]))
	elif temp_resname in other.resnames:
		other_index.append(int(line[0]))

nPro_res = len(pro_index)
nRNA_res = len(rna_index)
nOther_res = len(other_index)

res_list = np.zeros(nPro_res)
for i in range(nPro_res):
	res_list[i] = pro_index[i] + 168

bool_RNA = True
if nRNA_res == 0:
	bool_RNA = False

bool_other = True
if nOther_res == 0:
	bool_other = False

j = 0
nWindows = 0
for i in range(start,end,step):
	nWindows += 1

# LOOPING THROUGH ALL WINDOWS AND PLOTTING THE INDIVIDUAL DATA SETS
j = 0
window_count = 0
window_data = np.zeros((nSelections,nWindows),dtype=np.float64)
for i in range(start,end,step):
	j += step
	change_dir('%03d.%03d.RMSF' %(i,j))
	temp_data = np.loadtxt('%03d.%03d.%s.rmsf.dat' %(i,j,system))
	# STORING DATA OF EACH WINDOW TO BE USED LATER TO PLOT AVG AND STDEV OF ALL RMSF VALUES
	window_data[window_count] = temp_data

	plt.figure(1)
	plt.plot(res_list[:],temp_data[pro_index[0]:pro_index[-1]+1])

	plt.figure(2)
	plot_1d(res_list[:],temp_data[pro_index[0]:pro_index[-1]+1],'k','Protein Residue Number','RMSF','%03d.%03d.%s'%(i,j,system),'pro_rmsf',yunits='$\AA$',x_lim=(res_list[0],res_list[-1]),y_lim=(0,7),plt_title='RMSF of Protein Residues (Traj %03d - %03d)' %(i,j))

	if bool_RNA:
		plt.figure(3)
		plt.plot(range(nRNA_res),temp_data[rna_index[0]:rna_index[-1]+1])

		plt.figure(4)
		plot_1d(range(nRNA_res),temp_data[rna_index[0]:rna_index[-1]+1],'k','RNA Residue Selections','RMSF','%03d.%03d.%s'%(i,j,system),'rna_rmsf',marker_style='.',yunits='$\AA$',x_lim=(0,nRNA_res-1),y_lim=(0,7),plt_title='RMSF of RNA Residue Selections (Traj %03d - %03d)' %(i,j))

	if bool_other:
		plt.figure(5)
		plt.plot(range(nOther_res),temp_data[other_index[0]:])
		
		plt.figure(6)
		plot_1d(range(nOther_res),temp_data[other_index[0]:],'k','ATP Binding Substrate Selections','RMSF','%03d.%03d.%s'%(i,j,system),'atp_rmsf',marker_style='.',yunits='$\AA$',x_lim=(0,nOther_res-1),y_lim=(0,7),plt_title='RMSF of ATP Binding Pocket Substrate Selections (Traj %03d - %03d)' %(i,j))

	change_dir('..')

plt.figure(1)
plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
plt.xlabel('Protein Residue Number')
plt.ylabel('RMSF ($\AA$)')
plt.ylim((0,7))
plt.xlim((res_list[0],res_list[-1]))
plt.savefig('protein_rmsf_all_windows.png',dpi=300)
plt.close()

if bool_RNA:
	plt.figure(3)
	plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
	plt.xlabel('RNA Residue Selections',size=14)
	plt.ylabel('RMSF ($\AA$)',size=14)
	plt.xlim((0,nRNA_res-1))
	plt.ylim((0,7.0))
	plt.title('Window Comparison of RNA Residue RMSF',size=14)
	plt.savefig('rna_rmsf_all_windows.png',dpi=300)
	plt.close()

if bool_other:
	plt.figure(5)
	plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
	plt.xlabel('ATP Substrate Residue Selections',size=14)
	plt.ylabel('RMSF ($\AA$)',size=14)
	plt.xlim((0,nOther_res-1))
	plt.ylim((0,7.0))
	plt.title('Window Comparison of ATP Substrate Residue RMSF',size=14)
	plt.savefig('atp_sub_rmsf_all_windows.png',dpi=300)
	plt.close()

# AVERAGING AND CALCING THE STDEV OF ALL WINDOW'S RMSF VALUES
avg_data = np.mean(window_data,axis=1)
std_data = np.std(window_data,axis=1)

plt.errorbar(res_list[:],avg_data[pro_index[0]:pro_index[-1]+1],yerr=std_data[pro_index[0]:pro_index[-1]+1],marker='.')
plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
plt.xlabel('Protein Residue Number')
plt.ylabel('RMSF ($\AA$)')
plt.ylim((0,7))
plt.xlim((res_list[0],res_list[-1]))
plt.savefig('%03d.%03d.protein_rmsf_avg_errorbars.png'%(start,end),dpi=300)
plt.close()

if bool_RNA:
	plt.errorbar(range(nRNA_res),avg_data[rna_index[0]:rna_index[-1]+1],yerr=std_data[rna_index[0]:rna_index[-1]+1],marker='.')
	plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
	plt.xlabel('RNA Residue Selections')
	plt.ylabel('RMSF ($\AA$)')
	plt.ylim((0,7))
	plt.xlim((0,nRNA_res-1))
	plt.savefig('%03d.%03drna_rmsf_avg_errorbars.png'%(start,end),dpi=300)
	plt.close()

if bool_other:
	plt.errorbar(range(nOther_res),avg_data[other_index[0]:],yerr=std_data[other_index[0]:],marker='.')
	plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
	plt.xlabel('ATP Substrate Residue Selections')
	plt.ylabel('RMSF ($\AA$)')
	plt.ylim((0,7))
	plt.xlim((0,nOther_res-1))
	plt.savefig('%03d.%03d.atp_sub_rmsf_avg_errorbars.png'%(start,end),dpi=300)
	plt.close()

change_dir('021.150.RMSF')
data = np.loadtxt('021.150.%s.rmsf.dat' %(system))
plot_1d(res_list[:],data[pro_index[0]:pro_index[-1]+1],'k','Protein Residue Number','RMSF','021.150.%s'%(system),'pro_rmsf',yunits='$\AA$',x_lim=(res_list[0],res_list[-1]),y_lim=(0,7),plt_title='RMSF of Protein Residues (Traj 021 - 150)')

if nOther_res != 0:
	plot_1d(range(nOther_res),data[other_index[0]:],'k','Substrate Residue Selections','RMSF', '021.150.%s'%(system),'sub_rmsf',marker_style='.',yunits='$\AA$',x_lim=(0,nOther_res-1),y_lim=(0,4),plt_title='RMSF of Substrate Residue Selections (Traj 021 - 150)')

