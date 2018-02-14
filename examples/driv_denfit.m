shells = readmol('water.cc-pvdz.mol');
shells_aux = readmol('water.cc-pvdz-ri.mol');

nshell = length(shells);
naux   = length(shells_aux);

am = [shells.am];
nf = (am+1).*(am+2)/2; % number of functions
nftot = sum(nf);
off = [1 cumsum(nf)+1];

am_aux = [shells_aux.am];
nf_aux = (am_aux+1).*(am_aux+2)/2; % number of functions
nftot_aux = sum(nf_aux);
ofx = [1 cumsum(nf_aux)+1];

% unit shell is a single primitive with orbital exponent zero
unit_shell.atom_ind = -1;
unit_shell.atom_sym = 'X';
unit_shell.am       = 0;
unit_shell.nprim    = 1;
unit_shell.x        = 0;
unit_shell.y        = 0;
unit_shell.z        = 0;
unit_shell.alpha    = 0;
unit_shell.coef     = 1;

% calculate 3-center density fitting integrals
pqA = zeros(nftot,nftot,nftot_aux);
for i = 1:nshell
for j = 1:nshell
for k = 1:naux
  v = calculate_eri(shells(i), shells(j), shells_aux(k), unit_shell);
  v = reshape(v, nf_aux(k), nf(j), nf(i));
  v = permute(v, [3 2 1]);

  pqA(off(i):off(i+1)-1,...
      off(j):off(j+1)-1,...
      ofx(k):ofx(k+1)-1) = v;
end
end
end

% calculate the Coulomb metric matrix
J = zeros(nftot_aux,nftot_aux);
for i = 1:naux
for j = 1:naux
  v = calculate_eri(shells_aux(i), unit_shell, shells_aux(j), unit_shell);
  v = reshape(v, nf_aux(j), nf_aux(i));
  v = permute(v, [2 1]);
  J(ofx(i):ofx(i+1)-1, ofx(j):ofx(j+1)-1) = v;
end
end

% apply inverse square root
Jinvsqrt = inv(sqrtm(J));

% now check the answer...
% below tensor contractions can be much faster by unfolding into matrices...

% compute Temp = pqA*Jinvsqrt
Temp = zeros(size(pqA));
for i = 1:nftot
for j = 1:nftot
for k = 1:nftot_aux
  t = 0;
  for l = 1:nftot_aux
    t = t + pqA(i,j,l)*Jinvsqrt(l,k);
  end
  Temp(i,j,k) = t;
end
end
end

% compute Out = Temp*Temp'
Out = zeros(nftot,nftot,nftot,nftot);
for i = 1:nftot
for j = 1:nftot
for k = 1:nftot
for l = 1:nftot
  t = 0;
  for m = 1:nftot_aux
    t = t + Temp(i,j,m)*Temp(l,k,m);
  end
  Out(i,j,k,l) = t;
end
end
end
end

% calculate ERI for comparison
a = zeros(nftot,nftot,nftot,nftot);
for i = 1:nshell
for j = 1:nshell
for k = 1:nshell
for l = 1:nshell
  % not yet handling symmetries
  % THIS VERSION REVERSES SHELL ORDER AND DOES NOT NEED PERMUTATION
  v = calculate_eri(shells(l), shells(k), shells(j), shells(i));
  v = reshape(v, nf(i), nf(j), nf(k), nf(l));
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

norm(Out(:)-a(:))/norm(a(:))

