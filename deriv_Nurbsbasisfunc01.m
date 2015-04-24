function [R,dR] = deriv_Nurbsbasisfunc01(p,i,u,U,q,j,v,V,CP)
% returns the derivative of NURBS basis functions w.r.t u,v
% i,j are the knot span indices. if unknown, insert 0!
% dR is a vector of length (p+1)*(q+1)

if (i==0); i = findspan(u,U,length(CP(:,1,1)));  end
if (j==0); j = findspan(v,V,length(CP(1,:,1)));  end

ne = (p+1)*(q+1);      % Control Points per element
N = deriv(i,p,u,U);    % basisfunc in u
M = deriv(j,q,v,V);    % basisfunc in v

%  Form basis functions R and derivatives dR/du and dR/dv
k = 0;
sum = 0;
dsum = zeros(2,1);

for c = 0:q
  for b = 0:p
    k = k+1;
    % basis functions (only needed for the calc. of the derivatives)
    R(k) = N(1,b+1)*M(1,c+1)*CP(i-p+b,j-q+c,4);
    sum = sum + R(k);
            
    % derivatives
    dR(k,1) = N(2,b+1)*M(1,c+1)*CP(i-p+b,j-q+c,4);
    dsum(1)  = dsum(1) + dR(k,1);
    dR(k,2) = N(1,b+1)*M(2,c+1)*CP(i-p+b,j-q+c,4);
    dsum(2)  = dsum(2) + dR(k,2);
  end
end

    % divide through by sum
for k = 1:ne
    dR(k,1) = dR(k,1)/sum - (R(k)*dsum(1))/(sum^2);
    dR(k,2) = dR(k,2)/sum - (R(k)*dsum(2))/(sum^2);
  
    R(k) = R(k)/sum;
end