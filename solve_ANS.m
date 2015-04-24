function [d,fs] = solve_ANS(p,U,q,V,r,W,CP,E,nue,ngauss,f,rb)
% 1. calculate stiffness matrix
% 2. reduce stiffness matrix and load vector for fixed supports
% 3. calculate reduced displacement vector
% 4. assemble complete displacement vector d
% 5. calculate complete load vector fs
% For input parameter description see input file
% J. Kiendl

% 1. stiffness matrix
K = stiff_mat_ANS(p,U,q,V,r,W,CP,E,nue,ngauss);

% 2. reduce stiffness matrix and load vector
Kred = K;
fred = f;
for i = length(rb):-1:1
  Kred(:,rb(i)) = [];
  Kred(rb(i),:) = [];
  fred(rb(i)) = [];
end

% 3. reduced displacement vector
dred = Kred\fred;

% 4. assemble complete displacement vector
d = zeros(length(f),1);
i=1;
k=1;
for l = 1:length(f)
  if (i<=length(rb))
    if (l==rb(i))
      d(l)=0;
      i=i+1;
    else
      d(l)=dred(k);
      k=k+1;
    end
  else
    d(l)=dred(k);
    k=k+1;
  end
end

% 5. calculate complete load vector
fs=K*d;