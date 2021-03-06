%shells = readmol('water.sto-3g.mol');
%shells = readmol('mywater.mol');
%shells = readmol('water.cc-pvdz.mol');
%shells = readmol('graphene12.cc-pvdz.mol');
shells = readmol('alkane20.sto-3g.mol');

nshell = length(shells);

v = calculate_eri(shells(1), shells(2), shells(3), shells(4));

am = [shells.am];
nf = (am+1).*(am+2)/2; % number of functions
nftot = sum(nf);
off = [1 cumsum(nf)+1];

a = zeros(nftot,nftot,nftot,nftot);

for i = 1:nshell
for j = 1:nshell
for k = 1:nshell
for l = 1:nshell
  % not yet handling symmetries
  v = calculate_eri(shells(i), shells(j), shells(k), shells(l));
  %v = reshape(v, nf(i), nf(j), nf(k), nf(l));
  v = reshape(v, nf(l), nf(k), nf(j), nf(i));
  v = permute(v, [4 3 2 1]);
  a(off(i):off(i+1)-1,...
    off(j):off(j+1)-1,...
    off(k):off(k+1)-1,...
    off(l):off(l+1)-1) = ...
  a(off(i):off(i+1)-1,...
    off(j):off(j+1)-1,...
    off(k):off(k+1)-1,...
    off(l):off(l+1)-1) + v;
end
end
end
end

a0 = a;

