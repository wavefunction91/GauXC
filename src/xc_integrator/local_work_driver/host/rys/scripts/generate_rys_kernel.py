from itertools import product
import pyexpander.lib as expander
from math import ceil

s_targets = [ (0,0,0) ]
p_targets = [ (1,0,0), (0,1,0), (0,0,1) ]
d_targets = [ (2,0,0), (1,1,0), (1,0,1), (0,2,0), (0,1,1), (0,0,2) ]

sph_targets = [ s_targets, p_targets, d_targets ]

lA=1
lB=1

int_targets = []
for bra,ket in product( sph_targets[lA], sph_targets[lB] ):
    int_targets.append([bra,ket])

rys_targets  = []
uniq_targets = set() 
for bra,ket in int_targets:
    x_targets = (bra[0],ket[0])
    y_targets = (bra[1],ket[1])
    z_targets = (bra[2],ket[2])
    rys_targets.append([x_targets,y_targets,z_targets])
    uniq_targets.add(x_targets)
    uniq_targets.add(y_targets)
    uniq_targets.add(z_targets)
uniq_targets = list(uniq_targets)
uniq_targets.sort()


for u_target in uniq_targets:
    targets = [ target for target in rys_targets if u_target in target ]

nroot = ceil((lA+lB)/2.0) + 1
exp_dict = {'rys_targets':rys_targets, 'uniq_targets':uniq_targets, 'lA':lA, 'lB':lB, 'nroot':nroot }
expander.expandFile( 'rys_kernel_template.hpp', external_definitions=exp_dict, auto_indent=False )
