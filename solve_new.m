function d = solve_new(p,U,q,V,r,W,CP,E,nue,ngauss,f,rb)
% Instead of throwing rows and columns out of K, zero them out and put 1 on
% diagonal

% 1. stiffness matrix
K = stiff_mat(p,U,q,V,r,W,CP,E,nue,ngauss);

% 2. reduce stiffness matrix and load vector
for i = length(rb):-1:1
  K(:,rb(i)) = 0;
  K(rb(i),:) = 0;
  K(rb(i),rb(i)) = 1;
  f(rb(i)) = 0;
end

% 3. reduced displacement vector
d = K\f;
