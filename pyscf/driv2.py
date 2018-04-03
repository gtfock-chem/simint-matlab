import numpy as np
from pyscf import gto, scf

mol = gto.M(atom='''O  0.1009671920798554   0.07139585036771477 0.0
                    H  0.03308698470939767 -1.746407186150188   0.0,
                    H -1.6355095931157138   0.6133032301203728  0.0''',
            basis='sto3g')

# print(mol.bas_ctr_coeff(0))

ovl_mat = mol.intor('int1e_ovlp_cart')
print(ovl_mat)

kin_mat = mol.intor('int1e_kin_cart')
print(kin_mat)

nuc_mat = mol.intor('int1e_nuc_cart')
print(nuc_mat)

mf = scf.RHF(mol)
mf.kernel()

hcore = scf.hf.get_hcore(mol)
print(hcore)

guess = scf.hf.init_guess_by_atom(mol)
print(guess)

# SAD guesses
mol = gto.M(atom='''C  0 0 0''', basis='ccpvdz');
mol = gto.M(atom='''S  0 0 0''', basis='ccpvdz');
mol = gto.M(atom='''H  0 0 0,
                    H  1000 0 0''', basis='ccpvdz');
mol = gto.M(atom='''O  0 0 0''', basis='ccpvdz');
mf = scf.UHF(mol)
mf.kernel()
guess = scf.hf.init_guess_by_atom(mol)
print(guess.size);
print(guess)


