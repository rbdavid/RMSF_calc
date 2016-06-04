
sel = []
sel.append(['residue','protein and not name H*'])
sel.append(['backbone','backbone'])
sel.append(['sidechain','protein and not backbone and not (resname GLY or name H*)'])
sel.append(['C_alpha','protein and name CA'])

# CHECK THAT ATP, ADP, PI, and NUCLEIC RESIDUES ARE BEING ANALYZED... IF NOT CREATE ATOM SELECTIONS TO DO SO... 

