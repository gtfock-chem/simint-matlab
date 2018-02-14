% level 1 interface example

tau_screen = 1.e-11;

basis = gen_basis('../eri-matlab/bases/sto-3g.gbs', '../eri-matlab/molecules/water.xyz');
%basis = gen_basis('../eri-matlab/cc-pvdz.gbs', '../eri-matlab/molecules/water.xyz');

integral = gen_erd(basis);

num_atoms     = get_num_atoms(basis);
num_shells    = get_num_shells(basis);
num_functions = get_num_funcs(basis);

% what kind of indexing is used?
% dummy = get_shell_atom_idx(basis, 1); % what does it return?
% dummy = get_shell_start_idx(basis, 1);

% start of each shell
start = ones(num_shells+1,1);
for i = 1:num_shells
  start(i+1) = start(i) + get_shell_dim(basis, i);
end
if (start(num_shells+1)-1 ~= num_functions)
  error('error initializing start array');
end

% screening matrix
screenmat  = zeros(num_shells, num_shells);
screenmat2 = zeros(num_shells, num_shells);
for M = 1:num_shells
for N = 1:num_shells
  % could exploit symmetry (M,N,M,N) is same as (N,M,N,M)
  block = compute_shell_quartet(basis, integral, M, N, M, N);
  Mdim = get_shell_dim(basis, M);
  Ndim = get_shell_dim(basis, N);
  % compute shell pair value
  for i = 1:Mdim
    for j = 1:Ndim
      ele = block(i,j,i,j);
      if (ele < 0), error('whoa'); end
      screenmat(M,N) = max(screenmat(M,N), ele);
    end
  end
  screenmat2(M,N) = max(abs(block(:)));
end
end
if (norm(screenmat-screenmat2,'fro') ~= 0), error('screenmat error'); end

% ERI tensor
a = zeros(num_functions, num_functions, num_functions, num_functions);
for M = 1:num_shells
for N = 1:num_shells
for P = 1:num_shells
for Q = 1:num_shells
  % undone: exploit symmetries
  % Schwarz screening
  if (screenmat(M,N)*screenmat(P,Q) >= tau_screen*tau_screen)
    block = compute_shell_quartet(basis, integral, M, N, P, Q);
    a(start(M):start(M+1)-1,...
      start(N):start(N+1)-1,...
      start(P):start(P+1)-1,...
      start(Q):start(Q+1)-1) = ...
    a(start(M):start(M+1)-1,...
      start(N):start(N+1)-1,...
      start(P):start(P+1)-1,...
      start(Q):start(Q+1)-1) + block;
  end
end
end
end
end

destroy_basis(basis);
destroy_erd(integral);

