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
nRNA_res = 0
nOther_res = 0
nSelections = 0
for line in selections:
	line = line.split()
	temp_resname = line[1]
	nSelections += 1
	if temp_resname in pro.resnames:
		pro_index.append(int(line[0]))
	elif temp_resname in rna.resnames:
		nRNA_res += 1
	elif temp_resname in other.resnames:
		nOther_res += 1

nPro_res = len(pro_index)

pro_res_list = np.zeros(nPro_res)
for i in range(nPro_res):
	pro_res_list[i] = pro_index[i] 
#	pro_res_list[i] = pro_index[i] + 168

bool_RNA = False
if nRNA_res != 0:
	bool_RNA = True
	rna_min_index = nPro_res
	rna_max_index = rna_min_index + nRNA_res

bool_other = False
if nOther_res != 0:
	bool_other = True
	other_min_index = nPro_res+nRNA_res
	other_max_index = other_min_index + nOther_res

j = start-1
nWindows = 0
for i in range(start,end,step):
	nWindows += 1

# LOOPING THROUGH ALL WINDOWS AND PLOTTING THE INDIVIDUAL DATA SETS
j = start-1
window_count = 0
window_data = np.zeros((nSelections,nWindows),dtype=np.float64)
for i in range(start,end,step):
	j += step
	change_dir('%03d.%03d.RMSF' %(i,j))
	temp_data = np.loadtxt('%03d.%03d.%s.rmsf.dat' %(i,j,system))
	# STORING DATA OF EACH WINDOW TO BE USED LATER TO PLOT AVG AND STDEV OF ALL RMSF VALUES
	window_data[:,window_count] = temp_data

	plt.figure(1)
	plt.plot(pro_res_list[:],temp_data[pro_index[0]:pro_index[-1]+1])

	plt.figure(2)
	plot_1d(pro_res_list[:],temp_data[pro_index[0]:pro_index[-1]+1],'k','Protein Residue Number','RMSF','%03d.%03d.%s'%(i,j,system),'pro_rmsf',yunits='$\AA$',x_lim=(pro_res_list[0],pro_res_list[-1]),y_lim=(0,7),plt_title='RMSF of Protein Residues (Traj %03d - %03d)' %(i,j))

	if bool_RNA:
		plt.figure(3)
		plt.plot(range(nRNA_res),temp_data[rna_min_index:rna_max_index])

		plt.figure(4)
		plot_1d(range(nRNA_res),temp_data[rna_min_index:rna_max_index],'k','RNA Residue Selections','RMSF','%03d.%03d.%s'%(i,j,system),'rna_rmsf',marker_style='.',yunits='$\AA$',x_lim=(0,nRNA_res-1),y_lim=(0,4),plt_title='RMSF of RNA Residue Selections (Traj %03d - %03d)' %(i,j))

	if bool_other:
		plt.figure(5)
		plt.plot(range(nOther_res),temp_data[other_min_index:other_max_index])
		
		plt.figure(6)
		plot_1d(range(nOther_res),temp_data[other_min_index:other_max_index],'k','ATP Binding Substrate Selections','RMSF','%03d.%03d.%s'%(i,j,system),'atp_rmsf',marker_style='.',yunits='$\AA$',x_lim=(0,nOther_res-1),y_lim=(0,4),plt_title='RMSF of ATP Binding Pocket Substrate Selections (Traj %03d - %03d)' %(i,j))

	window_count += 1
	change_dir('..')

plt.figure(1)
plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
plt.xlabel('Protein Residue Number')
plt.ylabel('RMSF ($\AA$)')
plt.ylim((0,7))
plt.xlim((pro_res_list[0],pro_res_list[-1]))
plt.title('Window Comparison of Protein Residue RMSF - %s'%(system),size=14)
plt.savefig('%03d.%03d.%s.protein_rmsf_all_windows.png'%(start,end,system),dpi=300)
plt.close()

if bool_RNA:
	plt.figure(3)
	plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
	plt.xlabel('RNA Residue Selections',size=14)
	plt.ylabel('RMSF ($\AA$)',size=14)
	plt.xlim((0,nRNA_res-1))
	plt.ylim((0,4))
	plt.title('Window Comparison of RNA Residue RMSF - %s'%(system),size=14)
	plt.savefig('%03d.%03d.%s.rna_rmsf_all_windows.png'%(start,end,system),dpi=300)
	plt.close()

if bool_other:
	plt.figure(5)
	plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
	plt.xlabel('ATP Substrate Residue Selections',size=14)
	plt.ylabel('RMSF ($\AA$)',size=14)
	plt.xlim((0,nOther_res-1))
	plt.ylim((0,4))
	plt.title('Window Comparison of ATP Substrate Residue RMSF- %s'%(system),size=14)
	plt.savefig('%03d.%03d.%s.atp_sub_rmsf_all_windows.png'%(start,end,system),dpi=300)
	plt.close()

# AVERAGING AND CALCING THE STDEV OF ALL WINDOW'S RMSF VALUES
avg_data = np.mean(window_data,axis=1)
std_data = np.std(window_data,axis=1)

plt.errorbar(pro_res_list[:],avg_data[pro_index[0]:pro_index[-1]+1],yerr=std_data[pro_index[0]:pro_index[-1]+1],color='k',lw=2,elinewidth=1,ecolor='r')
plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
plt.xlabel('Protein Residue Number')
plt.ylabel('RMSF ($\AA$)')
plt.ylim((0,7))
plt.xlim((pro_res_list[0],pro_res_list[-1]))
plt.title('Avg w/ St Dev eror bars of Window RMSF of Protein Residues - %s'%(system))
plt.savefig('%03d.%03d.%s.protein_rmsf_avg_errorbars.png'%(start,end,system),dpi=300)
plt.close()

if bool_RNA:
	plt.errorbar(range(nRNA_res),avg_data[rna_min_index:rna_max_index],yerr=std_data[rna_min_index:rna_max_index],marker='.',color='k',lw=2,elinewidth=1,ecolor='r')
	plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
	plt.xlabel('RNA Residue Selections')
	plt.ylabel('RMSF ($\AA$)')
	plt.ylim((0,4))
	plt.xlim((0,nRNA_res-1))
	plt.title('Avg w/ St Dev eror bars of Window RMSF of RNA Selections - %s'%(system))
	plt.savefig('%03d.%03d.%s.rna_rmsf_avg_errorbars.png'%(start,end,system),dpi=300)
	plt.close()

if bool_other:
	plt.errorbar(range(nOther_res),avg_data[other_min_index:other_max_index],yerr=std_data[other_min_index:other_max_index],marker='.',color='k',lw=2,elinewidth=1,ecolor='r')
	plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
	plt.xlabel('ATP Substrate Residue Selections')
	plt.ylabel('RMSF ($\AA$)')
	plt.ylim((0,4))
	plt.xlim((0,nOther_res-1))
	plt.title('Avg w/ St Dev eror bars of Window RMSF of ATP Substrates - %s'%(system))
	plt.savefig('%03d.%03d.%s.atp_sub_rmsf_avg_errorbars.png'%(start,end,system),dpi=300)
	plt.close()

change_dir('021.150.RMSF')
data = np.loadtxt('021.150.%s.rmsf.dat' %(system))
plot_1d(pro_res_list[:],data[pro_index[0]:pro_index[-1]+1],'k','Protein Residue Number','RMSF','021.150.%s'%(system),'pro_rmsf',yunits='$\AA$',x_lim=(pro_res_list[0],pro_res_list[-1]),y_lim=(0,7),plt_title='RMSF of Protein Residues (Traj 021 - 150)')

if nRNA_res != 0:
	plot_1d(range(nRNA_res),data[rna_min_index:rna_max_index],'k','RNA Residue Selections','RMSF', '021.150.%s'%(system),'rna_rmsf',marker_style='.',yunits='$\AA$',x_lim=(0,nRNA_res-1),y_lim=(0,4),plt_title='RMSF of RNA Residue Selections (Traj 021 - 150)')

if nOther_res != 0:
	plot_1d(range(nOther_res),data[other_min_index:other_max_index],'k','ATP Substrate Selections','RMSF', '021.150.%s'%(system),'sub_rmsf',marker_style='.',yunits='$\AA$',x_lim=(0,nOther_res-1),y_lim=(0,4),plt_title='RMSF of ATP Substrate Selections (Traj 021 - 150)')

