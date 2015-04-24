function [d,fs] = solve_closed(p,U,q,V,r,W,CP,E,nue,ngauss,f,rb,ms)
% solution for a closed geometry;  ms = master-slaves-matrix
%
% 1. calculate stiffness matrix
% 2. reduce stiffness matrix and load vector for fixed supports
% 3. calculate reduced displacement vector
% 4. assemble complete displacement vector d
% 5. calculate complete load vector fs
% For input parameter description see input file

% 1. stiffness matrix
K = stiff_mat(p,U,q,V,r,W,CP,E,nue,ngauss);

% 2. reduce stiffness matrix and load vector
einsen=ones(1,length(K));  R=diag(einsen);
for i = length(ms(:,1)):-1:1
  R(:,ms(i,2)) = R(:,ms(i,2)) + R(:,ms(i,1));
  R(:,ms(i,1)) = [];
  for j = length(rb):-1:1    % adapt rb to ms
    if (rb(j) > ms(i,1));      rb(j) = rb(j)-1;
    elseif (rb(j) == ms(i,1)); rb(j) = [];  break;
    elseif (rb(j) < ms(i,1));  break;
    end
  end
end
Kred = K;
fred = f;
Kred = R'*Kred*R;
fred = R'*fred;
for i = length(rb):-1:1
  Kred(:,rb(i)) = [];
  Kred(rb(i),:) = [];
  fred(rb(i)) = [];
end

% 3. reduced displacement vector
dred = Kred\fred;

% 4. assemble complete displacement vector
d = zeros(length(f)-length(ms),1);
i=1;
k=1;
for l = 1:length(d)
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
d=R*d;
% 5. calculate complete load vector
fs=K*d;