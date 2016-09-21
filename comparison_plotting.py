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

u = MDAnalysis.Universe(pdb_file)
pro = u.select_atoms('protein')
other = u.select_atoms('nucleic or resname A5 U5 A3 atp adp PHX MG')

selections = open(selection_file,'r')

pro_index = []
other_index = []
for line in selections:
	line = line.split()
	temp_resname = line[1]
	if temp_resname in pro.resnames:
		pro_index.append(int(line[0]))
	elif temp_resname in other.resnames:
		other_index.append(int(line[0]))

nPro_res = len(pro_index)
nOther_res = len(other_index)

res_list = np.zeros(nPro_res)
for i in range(nPro_res):
	res_list[i] = pro_index[i] + 168

j = 0
legend_list = []
for i in range(start,end,step):
	j += step
	change_dir('%03d.%03d.RMSF' %(i,j))
	data = np.loadtxt('%03d.%03d.%s.rmsf.dat' %(i,j,system))
	
	legend_list.append('Traj %03d - %03d' %(i,j)) 
	
	plt.figure(1)
	plt.plot(res_list[:],data[pro_index[0]:pro_index[-1]+1])

	plt.figure(2)
	plot_1d(res_list[:],data[pro_index[0]:pro_index[-1]+1],'k','Protein Residue Number','RMSF','%03d.%03d.%s'%(i,j,system),'pro_rmsf',yunits='$\AA$',x_lim=(res_list[0],res_list[-1]),y_lim=(0,7),plt_title='RMSF of Protein Residues (Traj %03d - %03d)' %(i,j))
	
	if nOther_res != 0:
		plt.figure(3)
		plt.plot(range(nOther_res),data[other_index[0]:])
		
		plt.figure(4)
		plot_1d(range(nOther_res),data[other_index[0]:],'k','Substrate Residue Selections','RMSF','%03d.%03d.%s'%(i,j,system),'sub_rmsf',yunits='$\AA$',x_lim=(0,nOther_res-1),y_lim=(0,4),plt_title='RMSF of Substrate Residue Selections (Traj %03d - %03d)' %(i,j))

	change_dir('..')

plt.figure(1)
plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
#leg = plt.legend(legend_list, bbox_to_anchor=(-0.05, 1.03, 1.1, .100), fontsize='10', loc=3, ncol=4, mode="expand", borderaxespad=0., markerscale=3,numpoints=1)
plt.xlabel('Protein Residue Number')
plt.ylabel('RMSF ($\AA$)')
plt.ylim((0,7))
plt.xlim((res_list[0],res_list[-1]))
plt.savefig('protein_rmsf_all_windows.png',dpi=300)
plt.close()

if nOther_res != 0:
	plt.figure(3)
	plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
	plt.xlabel('Substrate Residue Selections',size=14)
	plt.ylabel('RMSF ($\AA$)',size=14)
	plt.xlim((0,nOther_res-1))
	plt.ylim((0,4.0))
	plt.title('Window Comparison of Substrate Residue RMSF',size=14)
	plt.savefig('substrate_rmsf_all_windows.png',dpi=300)
	plt.close()

change_dir('021.150.RMSF')
data = np.loadtxt('021.150.%s.rmsf.dat' %(system))
plot_1d(res_list[:],data[pro_index[0]:pro_index[-1]+1],'k','Protein Residue Number','RMSF','021.150.%s'%(i,j,system),'pro_rmsf',yunits='$\AA$',x_lim=(res_list[0],res_list[-1]),y_lim=(0,7),plt_title='RMSF of Protein Residues (Traj 021 - 150)' %(i,j))
plot_1d(range(nOther_res),data[other_index[0]:],'k','Substrate Residue Selections','RMSF', '021.150.%s'%(i,j,system),'sub_rmsf',yunits='$\AA$',x_lim=(0,nOther_res-1),y_lim=(0,4),plt_title='RMSF of Substrate Residue Selections (Traj 021 - 150)' %(i,j))

