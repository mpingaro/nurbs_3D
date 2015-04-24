function d = solve2(p,U,q,V,r,W,CP,E,nue,ngauss,f,rb)
% 1. calculate stiffness matrix
% put ones on diagonal
% J. Kiendl

% 1. stiffness matrix
K = stiff_mat(p,U,q,V,r,W,CP,E,nue,ngauss);

for i = length(rb):-1:1
  K(:,rb(i)) = 0;
  K(rb(i),:) = 0;
  K(rb(i),rb(i)) = 1;
  f(rb(i)) = 0;
end

% 3. reduced displacement vector
d = K\f;
