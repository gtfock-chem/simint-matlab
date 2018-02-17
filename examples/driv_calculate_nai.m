%shells = readmol('/home/edmond/now/simint-matlab/simint-v0.7/test/dat/water.sto-3g.mol');
shells = readmol('/home/echow/simint-generator/outdir/test/dat/water.sto-3g.mol');

% construct data structures for atoms in the molecule

num_atoms = max([shells.atom_ind]);

atomic_nums = zeros(num_atoms,1);
xcoords     = zeros(num_atoms,1);
ycoords     = zeros(num_atoms,1);
zcoords     = zeros(num_atoms,1);

for i = 1:length(shells)
  ind = shells(i).atom_ind;
  atomic_nums(ind) = atomic_num_map(shells(i).atom_sym);
  xcoords(ind)     = shells(i).x;
  ycoords(ind)     = shells(i).y;
  zcoords(ind)     = shells(i).z;
end

if min(atomic_nums) < 1
  error('Missing atom information for at least one atom index.');
end

% v = calculate_nai(atomic_nums, xcoords, ycoords, zcoords, shells(3), shells(3))

% calculate core Hamiltonian

am = [shells.am];
nf = (am+1).*(am+2)/2; % number of functions
nftot = sum(nf);
off = [1 cumsum(nf)+1];

Hcore = zeros(nftot,nftot);
nshell = length(shells);
for i = 1:nshell
for j = 1:nshell
% v = calculate_ovlpi(shells(i), shells(j));
  v = calculate_kei(shells(i), shells(j));
  v = v + calculate_nai(atomic_nums, xcoords, ycoords, zcoords, ...
                    shells(i), shells(j));
% v = reshape(v, nf(i), nf(j)); % check; H2O is not a good enough test
  v = reshape(v, nf(j), nf(i));
  v = v';
  Hcore(off(i):off(i+1)-1,off(j):off(j+1)-1) = ...
  Hcore(off(i):off(i+1)-1,off(j):off(j+1)-1) + v;
end
end
Hcore
