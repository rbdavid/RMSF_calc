
sel = []
sel.append(['residue','not name H*'])
sel.append(['backbone','backbone'])
#sel.append(['sidechain_test','resname GLY and name ...' NOT GOING TO FUCKING WORK...)]

#sel.append(['sidechain','not backbone and not (resname GLY or name H*)'])
sel.append(['C_alpha','protein and name CA'])

# CHECK THAT ATP, ADP, PI, and NUCLEIC RESIDUES ARE BEING ANALYZED... IF NOT CREATE ATOM SELECTIONS TO DO SO... 

