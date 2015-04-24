function I = integral_area(p,q,r,i,j,k,u,v,w,U,V,W,CP,uvw,proj)
% evaluates R(u,v,w)*dA for an integration point
% multiply with gauss weights and map in the calling function
% Parameters:
%     i,j,k:     element indices
%     p,q,r:     polynomial degrees
%     u,v,w:     integration point coordinates
%     U,V,W:     knot vectors
%     uvw:       directions to integrate: 1=vw, 2=wu, 3=uv
%     proj:      if 0: integrate dA
%                if 1: integrate the projections of dA to the x,y,z plane

ne = (p+1)*(q+1)*(r+1);  % Control Points per element
Nu = deriv(i,p,u,U);     % basisfunc in u
Nv = deriv(j,q,v,V);     % basisfunc in v
Nw = deriv(k,r,w,W);     % basisfunc in w

%  Form basis functions R and derivatives dR/du and dR/dv
l = 0;
sum = 0;
dsum = zeros(3,1);

R = zeros(1,(p+1)*(q+1)*(r+1));
dRl = zeros((p+1)*(q+1)*(r+1),3);
for c = 0:r
  for b = 0:q
    for a = 0:p
      l = l+1;
      % basis functions
      R(l) = Nu(1,a+1)*Nv(1,b+1)*Nw(1,c+1)*CP(i-p+a,j-q+b,k-r+c,4);
      sum = sum + R(l);
            
      % derivatives
      dRl(l,1) = Nu(2,a+1)*Nv(1,b+1)*Nw(1,c+1)*CP(i-p+a,j-q+b,k-r+c,4);
      dsum(1)  = dsum(1) + dRl(l,1);
      dRl(l,2) = Nu(1,a+1)*Nv(2,b+1)*Nw(1,c+1)*CP(i-p+a,j-q+b,k-r+c,4);
      dsum(2)  = dsum(2) + dRl(l,2);
      dRl(l,3) = Nu(1,a+1)*Nv(1,b+1)*Nw(2,c+1)*CP(i-p+a,j-q+b,k-r+c,4);
      dsum(3)  = dsum(3) + dRl(l,3);
    end
  end
end

    % divide dRl through by sum
for l = 1:ne
    dRl(l,1) = dRl(l,1)/sum - (R(l)*dsum(1))/(sum^2);
    dRl(l,2) = dRl(l,2)/sum - (R(l)*dsum(2))/(sum^2);
    dRl(l,3) = dRl(l,3)/sum - (R(l)*dsum(3))/(sum^2);
    R(l) = R(l)/sum;
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

% dA is the norm of the cross product of the element's edges ds1, ds2
s1 = uvw+1;  if(s1>3); s1=s1-3; end;    ds1(1:3,1)=J(s1,:);
s2 = uvw+2;  if(s2>3); s2=s2-3; end;    ds2(1:3,1)=J(s2,:);
if     (proj == 0)
  dA = sqrt(cross(ds1,ds2)' * cross(ds1,ds2));
  I = zeros(p+1,q+1,r+1);  l = 0;
  for c = 0:r
    for b = 0:q
      for a = 0:p
        l = l + 1;
        I(a+1,b+1,c+1) = R(l)*dA;
      end
    end
  end
elseif (proj == 1)
% for the projections of dA onto the 3 planes reduce ds1,ds2 to the planes.
% The cross product has only one entry ~=0 which therefore is equal to dA(e)
  for e=1:3
    ds1_prj = ds1;    ds1_prj(e) = 0;
    ds2_prj = ds2;    ds2_prj(e) = 0;
    cross_prod = cross(ds1_prj,ds2_prj);
    dA(e) = cross_prod(e);
  end
  I = zeros(p+1,q+1,r+1,3);  l = 0;
  for c = 0:r
    for b = 0:q
      for a = 0:p
        l = l + 1;
        for e=1:3
          I(a+1,b+1,c+1,e) = R(l)*dA(e);
        end
      end
    end
  end
end