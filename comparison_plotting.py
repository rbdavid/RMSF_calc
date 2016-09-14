#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python

import numpy as np
import sys
import os 
from sel_list import *
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

sel = sys.argv[1]
start = int(sys.argv[2])
end = int(sys.argv[3])


#file_list = []
#file_list.append(['Apo','../../../AMBER_apo/RMSF/51_to_100/Apo.rmsf.%s.dat' %(sel),'steelblue'])
#file_list.append(['ATP','../../../AMBER_atp/RMSF/51_to_100/ATP.rmsf.%s.dat' %(sel),'cadetblue'])
#file_list.append(['ssRNA','../../../AMBER_ssrna/RMSF/51_to_100/ssRNA.rmsf.%s.dat' %(sel),'turquoise'])
#file_list.append(['ssRNA+ATP','../../../AMBER_ssrna_atp/RMSF/51_to_100/ssRNA+ATP.rmsf.%s.dat' %(sel),'forestgreen'])
#file_list.append(['ssRNA+ADP+Pi','../../../AMBER_ssrna_adp_pi/RMSF/51_to_100/ssRNA+ADP+Pi.rmsf.%s.dat' %(sel),'limegreen'])
#file_list.append(['ssRNA+ADP','../../../AMBER_ssrna_adp/RMSF/51_to_100/ssRNA+ADP.rmsf.%s.dat' %(sel),'orangered'])
#file_list.append(['ssRNA+Pi','../../../AMBER_ssrna_pi/RMSF/51_to_100/ssRNA+Pi.rmsf.%s.dat' %(sel),'crimson'])

nSys = len(file_list)
#legend_list = []
#for i in range(nSys):
#	legend_list.append(file_list[i][0])

flush = sys.stdout.flush

def ffprint(string):
        print '%s' %(string)
        flush()

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


