function B = B_matrix_tl(i,p,u,U,j,q,v,V,k,r,w,W,CP,d)
% Evaluates the stiffness matrix of an element for one integration point
% Input: p,q,r:        polynomial degrees
%        U,V,W:        knot vectors
%        CP:           control points
%        i,j,k:        element index in U,V,W
%        u,v,w:        NURBS coordinates of integration point
%        D:            material matrix
%        gwu,gwv,gww:  gaussian weights
% J. Kiendl

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


% strain
BL = zeros(6,3*ne);
for zz = 1:ne
  
  BL(1,3*(zz-1)+1)= dRd(zz,1);
  BL(1,3*(zz-1)+2)= 0.0;
  BL(1,3*(zz-1)+3)= 0.0;
  
  BL(2,3*(zz-1)+1)= 0.0;
  BL(2,3*(zz-1)+2)= dRd(zz,2);
  BL(2,3*(zz-1)+3)= 0.0;
  
  BL(3,3*(zz-1)+1)= 0.0;
  BL(3,3*(zz-1)+2)= 0.0;
  BL(3,3*(zz-1)+3)= dRd(zz,3);
  
  BL(4,3*(zz-1)+1)= dRd(zz,2);
  BL(4,3*(zz-1)+2)= dRd(zz,1);
  BL(4,3*(zz-1)+3)= 0.0;
  
  BL(5,3*(zz-1)+1)= 0.0;
  BL(5,3*(zz-1)+2)= dRd(zz,3);
  BL(5,3*(zz-1)+3)= dRd(zz,2);
  
  BL(6,3*(zz-1)+1)= dRd(zz,3);
  BL(6,3*(zz-1)+2)= 0.0;
  BL(6,3*(zz-1)+3)= dRd(zz,1);

end

L = zeros(3,3);
for nl = 1:ne
    L(1,1) = L(1,1) + dRd(nl,1)*d(3*(nl-1)+1);
    L(1,2) = L(1,2) + dRd(nl,2)*d(3*(nl-1)+1);
    L(1,3) = L(1,3) + dRd(nl,3)*d(3*(nl-1)+1);
    
    L(2,1) = L(2,1) + dRd(nl,1)*d(3*(nl-1)+2);
    L(2,2) = L(2,2) + dRd(nl,2)*d(3*(nl-1)+2);
    L(2,3) = L(2,3) + dRd(nl,3)*d(3*(nl-1)+2);

    L(3,1) = L(3,1) + dRd(nl,1)*d(3*(nl-1)+3);
    L(3,2) = L(3,2) + dRd(nl,2)*d(3*(nl-1)+3);
    L(3,3) = L(3,3) + dRd(nl,3)*d(3*(nl-1)+3);
end

BLL = zeros(6,3*ne);
for nl = 1:ne
  
  BLL(1,3*(nl-1)+1)= L(1,1)*dRd(nl,1);
  BLL(1,3*(nl-1)+2)= L(2,1)*dRd(nl,1);
  BLL(1,3*(nl-1)+3)= L(3,1)*dRd(nl,1);
  
  BLL(2,3*(nl-1)+1)= L(1,2)*dRd(nl,2);
  BLL(2,3*(nl-1)+2)= L(2,2)*dRd(nl,2);
  BLL(2,3*(nl-1)+3)= L(3,2)*dRd(nl,2);
  
  BLL(3,3*(nl-1)+1)= L(1,3)*dRd(nl,3);
  BLL(3,3*(nl-1)+2)= L(2,3)*dRd(nl,3);
  BLL(3,3*(zz-1)+3)= L(3,3)*dRd(nl,3);
  
  BLL(4,3*(nl-1)+1)= L(1,1)*dRd(nl,2)+L(1,2)*dRd(nl,1);
  BLL(4,3*(nl-1)+2)= L(2,1)*dRd(nl,2)+L(2,2)*dRd(nl,1);
  BLL(4,3*(nl-1)+3)= L(3,1)*dRd(nl,2)+L(3,2)*dRd(nl,1);
  
  BLL(5,3*(nl-1)+1)= L(1,2)*dRd(nl,3)+L(1,3)*dRd(nl,2);
  BLL(5,3*(nl-1)+2)= L(2,2)*dRd(nl,3)+L(2,3)*dRd(nl,2);
  BLL(5,3*(nl-1)+3)= L(3,2)*dRd(nl,3)+L(3,3)*dRd(nl,2);
  
  BLL(6,3*(nl-1)+1)= L(1,1)*dRd(nl,3)+L(1,3)*dRd(nl,1);
  BLL(6,3*(nl-1)+2)= L(2,1)*dRd(nl,3)+L(2,3)*dRd(nl,1);
  BLL(6,3*(nl-1)+3)= L(3,1)*dRd(nl,3)+L(3,3)*dRd(nl,1);

end

B = BL+BLL;

return