A = load('ERI');
A = reshape(A, 7,7,7,7);
S = load('S_ovlp');
S = reshape(S, 7,7);
H = load('H_core');
H = reshape(H, 7,7);
G = load('D_guess'); % not currently used (but may be needed if SCF does not converge)
G = reshape(G, 7,7);
nuc_ene = load('nuclear_energy');
nelectron = load('nelectron');
if round(nelectron/2) ~= nelectron/2
  error('nelectron is not even');
end
nocc = nelectron/2; % number of occupied orbitals (restricted Hartree-Fock)

% use symmetric orthogonalization, X = S^{-1/2}
[u d] = eig(S);
X = u*inv(sqrt(d))*u';

% set initial F to core Hamiltonian
F = H;

for iter = 1:20

  F = X'*F*X;

  % only need eigenvectors corresponding to occupied orbitals
  [C E] = eig(F);
  [dummy ind] = sort(diag(E));
  Chat = C(:,ind(1:nocc));
  Chat = X*Chat;

  % in practice, do not need to form D, but store the Chat factor only
  D = Chat*Chat';

  F = H;
  for mu=1:7
  for nu=1:7
    F(mu,nu) = F(mu,nu) ...
        + sum(sum(D .* (2*squeeze(A(mu,nu,:,:)) - squeeze(A(mu,:,:,nu)))));
  end
  end

  % calculate energy (may not need at every step)
  ene = sum(sum(D .* (H+F))) + nuc_ene;
  fprintf('%3d  %.15f\n', iter, ene);
end

% sample output for water in sto3g
% ...
%  1  -73.196954311244767
% ...
% 20  -74.945021299866966

