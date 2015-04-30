function epsilon = strain_tl(p,q,r,i,j,k,u,v,w,U,V,W,CP,d)
% Calculates B-Matrix as in stiff_mat_el. Then strain = B*d
% Input: p,q,r:     polynomial degrees
%        U,V,W:     knot vectors
%        CP:        control points
%        u,v,w:     NURBS coordinates sampling point
%        d:         element displacement vector

ne = (p+1)*(q+1)*(r+1);      % Control Points per element
N = deriv(i,p,u,U);    % basisfunc in u
M = deriv(j,q,v,V);    % basisfunc in v
O = deriv(k,r,w,W);    % basisfunc in w

%  Form basis functions R and derivatives dR/du, dR/dv and dR/dw
l = 0;
sum = 0;
dsum = zeros(3,1);
R = zeros(1,(p+1)*(q+1)*(r+1));
dRl = zeros((p+1)*(q+1)*(r+1),3);
for c = 0:r
  for b = 0:q
    for a = 0:p
      l = l+1;
      % basis functions (only needed for the calc. of the derivatives)
      R(l) = N(1,a+1)*M(1,b+1)*O(1,c+1)*CP(i-p+a,j-q+b,k-r+c,4);
      sum = sum + R(l);
            
      % derivatives
      dRl(l,1) = N(2,a+1)*M(1,b+1)*O(1,c+1)*CP(i-p+a,j-q+b,k-r+c,4);
      dsum(1)  = dsum(1) + dRl(l,1);
      dRl(l,2) = N(1,a+1)*M(2,b+1)*O(1,c+1)*CP(i-p+a,j-q+b,k-r+c,4);
      dsum(2)  = dsum(2) + dRl(l,2);
      dRl(l,3) = N(1,a+1)*M(1,b+1)*O(2,c+1)*CP(i-p+a,j-q+b,k-r+c,4);
      dsum(3)  = dsum(3) + dRl(l,3);
    end
  end
end

    % divide through by sum
for l = 1:ne
    dRl(l,1) = dRl(l,1)/sum - (R(l)*dsum(1))/(sum^2);
    dRl(l,2) = dRl(l,2)/sum - (R(l)*dsum(2))/(sum^2);
    dRl(l,3) = dRl(l,3)/sum - (R(l)*dsum(3))/(sum^2);
end

% Jacobi
J = zeros(3,3);
l = 0;
for c = 0:r
  for b = 0:q
    for a = 0:p
      l = l + 1;
      J(1,1) = J(1,1) + CP(i-p+a,j-q+b,k-r+c,1)*dRl(l,1);
      J(1,2) = J(1,2) + CP(i-p+a,j-q+b,k-r+c,2)*dRl(l,1);
      J(1,3) = J(1,3) + CP(i-p+a,j-q+b,k-r+c,3)*dRl(l,1);
      J(2,1) = J(2,1) + CP(i-p+a,j-q+b,k-r+c,1)*dRl(l,2);
      J(2,2) = J(2,2) + CP(i-p+a,j-q+b,k-r+c,2)*dRl(l,2);
      J(2,3) = J(2,3) + CP(i-p+a,j-q+b,k-r+c,3)*dRl(l,2);
      J(3,1) = J(3,1) + CP(i-p+a,j-q+b,k-r+c,1)*dRl(l,3);
      J(3,2) = J(3,2) + CP(i-p+a,j-q+b,k-r+c,2)*dRl(l,3);
      J(3,3) = J(3,3) + CP(i-p+a,j-q+b,k-r+c,3)*dRl(l,3);
    end
  end
end
dxidx = inv(J);
dRd = dRl*dxidx';

% Green-Lagrange strain
B = zeros(6,3*ne);
for zz = 1:ne
  
  B(1,3*(zz-1)+1)= dRd(zz,1)+0.5*dRd(zz,1)^2;
  B(1,3*(zz-1)+2)= 0.5*dRd(zz,1)^2;
  B(1,3*(zz-1)+3)= 0.5*dRd(zz,1)^2;
  
  B(2,3*(zz-1)+1)= 0.5*dRd(zz,2)^2;
  B(2,3*(zz-1)+2)= dRd(zz,2)+0.5*dRd(zz,2)^2;
  B(2,3*(zz-1)+3)= 0.5*dRd(zz,2)^2;
  
  B(3,3*(zz-1)+1)= 0.5*dRd(zz,3)^2;
  B(3,3*(zz-1)+2)= 0.5*dRd(zz,3)^2;
  B(3,3*(zz-1)+3)= dRd(zz,3)+0.5*dRd(zz,3)^2;
  
  B(4,3*(zz-1)+1)= dRd(zz,2)+dRd(zz,1)*dRd(zz,2);
  B(4,3*(zz-1)+2)= dRd(zz,1)+dRd(zz,1)*dRd(zz,2);
  B(4,3*(zz-1)+3)= dRd(zz,1)*dRd(zz,2);
  
  B(5,3*(zz-1)+1)= dRd(zz,2)*dRd(zz,3);
  B(5,3*(zz-1)+2)= dRd(zz,3)+dRd(zz,2)*dRd(zz,3);
  B(5,3*(zz-1)+3)= dRd(zz,2)+dRd(zz,2)*dRd(zz,3);
  
  B(6,3*(zz-1)+1)= dRd(zz,3)+dRd(zz,1)*dRd(zz,3);
  B(6,3*(zz-1)+2)= dRd(zz,1)*dRd(zz,3);
  B(6,3*(zz-1)+3)= dRd(zz,1)+dRd(zz,1)*dRd(zz,3);

end

% strain
epsilon = B*d;

return