#box length CHECK BASE.TXT!!
set box 100
#set filename sd14.mol2
#mol addfile $filename
#set filename pol200_e65_q10.pdb
#mol addfile $filename
display projection orthographic
display resetview
display depthcue off
color Display Background white
#B monomer
mol modselect 0 0 type  BB A2B A5B A3B A6B A4B
mol modcolor 0 0 ColorID 0
mol modstyle 0 0 VDW 0.5 100.0
mol addrep 0
#C monomer
mol modselect 1 0 type  BC A2C A3C A5C A4C A6C 
mol modcolor 1 0 ColorID 1
mol modstyle 1 0 VDW 0.5 100.0
mol addrep 0
#A1B bead is shared by both monomers
mol modselect 2 0 type A1B
mol modcolor 2 0 ColorID 9
mol modstyle 2 0 VDW 0.5 100.0
mol addrep 0
#A-A dimers
mol modselect 3 0 type  BA A1A A2A A5A A3A A6A A4A
mol modcolor 3 0 ColorID 19
mol modstyle 3 0 VDW 0.5 100.0
mol addrep 0
mol modselect 4 0 type TB
mol modcolor 4 0 ColorID 0
mol modstyle 4 0 VDW 0.7 100.0
mol addrep 0
mol modselect 5 0 type TC
mol modcolor 5 0 ColorID 1
mol modstyle 5 0 VDW 0.7 100.0
mol addrep 0
mol modselect 6 0 type TA
mol modcolor 6 0 ColorID 19
mol modstyle 6 0 VDW 0.7 100.0
mol addrep 0
#exluders
mol modselect 7 0 type X1 C0 C1 C2
mol modcolor 7 0 ColorID 2
mol modstyle 7 0 VDW 0.5 100.0
mol addrep 0
#polymer
mol modselect 8 0 type P
mol modcolor 8 0 ColorID 16
mol modstyle 8 0 VDW 0.4 100.0
mol addrep 0
#T-1 dimers
mol modselect 9 0 type B A1 A2 A3 A4 A5 A6
mol modcolor 9 0 ColorID 10
mol modstyle 9 0 VDW 0.5 100.0
mol addrep 0
mol modselect 10 0 type T 
mol modcolor 10 0 ColorID 10
mol modstyle 10 0 VDW 0.7 100.0
mol addrep 0
#set filename sd9.mol2
#mol addfile $filename

#periodic boundaries for mol2 files
molinfo top set {a b c alpha beta gamma} {$box $box $box 90 90 90}
menu main on
menu graphics on
animate goto end
animate speed 0.7
animate style once
