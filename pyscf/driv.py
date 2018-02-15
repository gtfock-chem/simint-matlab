import numpy as np
from pyscf import gto, scf

mol = gto.M(atom='''O 0.000000  0.000000 -0.061665  
                    H 0.000000 -0.711621  0.489331  
                    H 0.000000  0.711621  0.489331''',
            basis='sto3g')
mf = scf.RHF(mol)
mf.kernel()

f = open('S_ovlp', 'w')
mol.intor('int1e_ovlp').tofile(f, '\n')
f.close()

f = open('ERI', 'w')
mol.intor('int2e').tofile(f, '\n')
f.close()

f = open('H_core', 'w')
scf.hf.get_hcore(mol).tofile(f, '\n')
f.close()

f = open('D_guess', 'w')
scf.hf.init_guess_by_atom(mol).tofile(f, '\n')
f.close()

f = open('nuclear_energy', 'w')
mol.energy_nuc().tofile(f, '\n')
f.close()

f = open('nelectron', 'w')
f.write(str(mol.nelectron))
f.close()
