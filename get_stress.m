function S = get_stress(p,q,r,u,v,w,U,V,W,CP,d,E,nue)
% returns all stress components at a point u,v,w
% components 1-6: normal and shear stresses in x,y,z
%            7-9: principal stresses
%            10:  von Mises stress

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
epsilon = strain(p,q,r,i,j,k,u,v,w,U,V,W,CP,d_el);

% material matrix
D = E/((1+nue)*(1-2*nue))*[1-nue nue nue 0 0 0; nue 1-nue nue 0 0 0; nue nue 1-nue 0 0 0
                           0 0 0 (1-2*nue)/2 0 0; 0 0 0 0 (1-2*nue)/2 0; 0 0 0 0 0 (1-2*nue)/2];

sig = D*epsilon;

principal stresses from eigenvalues of stress tensor
S = [sig(1) sig(4) sig(6)
     sig(4) sig(2) sig(5)
     sig(6) sig(5) sig(3)];
 
s = eig(S);
s = sort(s,'descend');
sig(7:9) = s(1:3);
% von Mises stress
sig(10) = sqrt(sig(1)^2+sig(2)^2+sig(3)^2-sig(1)*sig(2)-sig(2)*sig(3)-sig(3)*sig(1)+3*sig(4)^2+3*sig(5)^2+3*sig(6)^2);