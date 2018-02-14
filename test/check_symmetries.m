function check_symmetries(a)
n = length(a);

d = zeros(1,8);
for i = 1:n
for j = 1:n
for k = 1:n
for l = 1:n
  d(1) = a(i,j,k,l) - a(i,j,k,l);
  d(2) = a(i,j,k,l) - a(i,j,l,k);
  d(3) = a(i,j,k,l) - a(j,i,k,l);
  d(4) = a(i,j,k,l) - a(j,i,l,k);
  d(5) = a(i,j,k,l) - a(k,l,i,j);
  d(6) = a(i,j,k,l) - a(k,l,j,i);
  d(7) = a(i,j,k,l) - a(l,k,i,j);
  d(8) = a(i,j,k,l) - a(l,k,j,i);
  if norm(d) > 1e-14
    disp(norm(d))
    d
  end
end
end
end
end

