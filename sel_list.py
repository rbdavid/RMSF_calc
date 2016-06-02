
sel_list = []
sel_list.append(['residue','protein and not name H*'])
sel_list.append(['backbone','backbone'])
sel_list.append(['sidechain','protein and not backbone and not name H*'])
sel_list.append(['C_alpha','protein and name CA'])

# CHECK THAT ATP, ADP, PI, and NUCLEIC RESIDUES ARE BEING ANALYZED... IF NOT CREATE ATOM SELECTIONS TO DO SO... 

