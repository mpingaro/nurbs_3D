function eps = get_strain(p,q,r,u,v,w,U,V,W,CP,d)
% returns all strain components at a point u,v,w
% components 1-6: normal and shear strains in x,y,z
%            7-9: principal strains

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

eps = strain(p,q,r,i,j,k,u,v,w,U,V,W,CP,d_el);

% principal strains from eigenvalues of strain tensor
E = [    eps(1) 0.5*eps(4) 0.5*eps(6)
     0.5*eps(4)     eps(2) 0.5*eps(5)
     0.5*eps(6) 0.5*eps(5)     eps(3)];
 
e = eig(E);
e = sort(e,'descend');
eps(7:9) = e(1:3);