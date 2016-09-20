#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:

# ----------------------------------------
# PREAMBLE:

from plotting_functions import *
import sys
import os 

start = int(sys.argv[1])
end = int(sys.argv[2])
step = int(sys.argv[3])

change_dir = os.chdir

for i in range(start,end,step):
	j += step
	change_dir('%03d.%03d.RMSF' %(i,j))


for i in range(nSys):
	data = np.loadtxt(file_list[i][1])
	nRes = len(data)
	res_list = np.zeros(nRes)
	for j in range(nRes):
		res_list[j] = j + 168
	plt.plot(res_list[:],data[:,1], c=file_list[i][2])

plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
leg = plt.legend(legend_list, bbox_to_anchor=(-0.05, 1.03, 1.1, .100), fontsize='10', loc=3, ncol=4, mode="expand", borderaxespad=0., markerscale=3,numpoints=1)
plt.xlabel('Residue Number')
plt.ylabel('RMSF ($\AA$)')
plt.ylim((0,4))
plt.xlim((res_list[0],res_list[-1]))

plt.savefig('%s.comparison.png' %(sel))
plt.close()


