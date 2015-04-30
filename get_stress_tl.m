function S = get_stress_tl(p,q,r,u,v,w,U,V,W,CP,d,E,nue)
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

epsilon = strain_tl(p,q,r,i,j,k,u,v,w,U,V,W,CP,d_el);

% material matrix
D = E/((1+nue)*(1-2*nue))*[1-nue nue nue 0 0 0; nue 1-nue nue 0 0 0; nue nue 1-nue 0 0 0
                           0 0 0 (1-2*nue)/2 0 0; 0 0 0 0 (1-2*nue)/2 0; 0 0 0 0 0 (1-2*nue)/2];

S = D*epsilon;

return