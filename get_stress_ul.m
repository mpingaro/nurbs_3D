function sig_v = get_stress_ul(p,q,r,u,v,w,U,V,W,CP,CP0,d,E,nue)
% returns all stress components at a point u,v,w

nu = length(CP(:,1,1,1));
nv = length(CP(1,:,1,1));
nw = length(CP(1,1,:,1));

i = findspan(u,U,nu);
j = findspan(v,V,nv);
k = findspan(w,W,nw);

% element displacement vector
d_el = zeros(3*(p+1)*(q+1)*(r+1),1);
l=1;
for c = k-r-1:k-1
  for b = j-q-1:j-1 
    for a = i-p-1:i-1
      d_el(l,1)   = d(3*(c*nu*nv+b*nu+a)+1);
      d_el(l+1,1) = d(3*(c*nu*nv+b*nu+a)+2);
      d_el(l+2,1) = d(3*(c*nu*nv+b*nu+a)+3);
      l=l+3;
    end
  end
end

F0 = jacobian(i,p,u,U,j,q,v,V,k,r,w,W,CP0);
Ft = jacobian(i,p,u,U,j,q,v,V,k,r,w,W,CP);
%F = Ft*inv(F0);
F =  Ft*(F0\blkdiag(1,1,1));
J = det(F);

epsilon = strain_ul(p,q,r,i,j,k,u,v,w,U,V,W,CP0,d_el);
%epsilon = strain(p,q,r,i,j,k,u,v,w,U,V,W,CP0,d_el);

% material matrix
D = E/((1+nue)*(1-2*nue))*[1-nue nue nue 0 0 0; nue 1-nue nue 0 0 0; nue nue 1-nue 0 0 0
                           0 0 0 (1-2*nue)/2 0 0; 0 0 0 0 (1-2*nue)/2 0; 0 0 0 0 0 (1-2*nue)/2];

S = D*epsilon;

St  = [S(1), S(4), S(6); S(4), S(2), S(5); S(6), S(5), S(3)];
sig = F*St*F'/J;

sig_v = [sig(1,1);sig(2,2);sig(3,3);sig(1,2);sig(2,3);sig(3,1)];

return