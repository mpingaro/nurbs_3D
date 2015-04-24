function BNL = B_NL_matrix(i,p,u,U,j,q,v,V,k,r,w,W,CP)

ne = (p+1)*(q+1)*(r+1);      % Control Points per element
N = deriv(i,p,u,U);          % basisfunc in u
M = deriv(j,q,v,V);          % basisfunc in v
O = deriv(k,r,w,W);          % basisfunc in w

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
dxdxi = zeros(3,3);
l = 0;
for c = 0:r
  for b = 0:q
    for a = 0:p
      l = l + 1;
      dxdxi(1,1) = dxdxi(1,1) + CP(i-p+a,j-q+b,k-r+c,1)*dRl(l,1);
      dxdxi(2,1) = dxdxi(2,1) + CP(i-p+a,j-q+b,k-r+c,2)*dRl(l,1);
      dxdxi(3,1) = dxdxi(3,1) + CP(i-p+a,j-q+b,k-r+c,3)*dRl(l,1);
      dxdxi(1,2) = dxdxi(1,2) + CP(i-p+a,j-q+b,k-r+c,1)*dRl(l,2);
      dxdxi(2,2) = dxdxi(2,2) + CP(i-p+a,j-q+b,k-r+c,2)*dRl(l,2);
      dxdxi(3,2) = dxdxi(3,2) + CP(i-p+a,j-q+b,k-r+c,3)*dRl(l,2);
      dxdxi(1,3) = dxdxi(1,3) + CP(i-p+a,j-q+b,k-r+c,1)*dRl(l,3);
      dxdxi(2,3) = dxdxi(2,3) + CP(i-p+a,j-q+b,k-r+c,2)*dRl(l,3);
      dxdxi(3,3) = dxdxi(3,3) + CP(i-p+a,j-q+b,k-r+c,3)*dRl(l,3);
    end
  end
end
dxidx = inv(dxdxi);
dRd = dRl*dxidx';

B = zeros(3,3*ne);
for i = 1:ne
  B(1,3*(i-1)+1)=dRd(i,1);
  B(2,3*(i-1)+1)=dRd(i,2);
  B(3,3*(i-1)+1)=dRd(i,3);
end
BNL = [B; zeros(3,1),B(:,1:end-1); zeros(3,2),B(:,1:end-2)];

return