function create_el_edges(p,q,U,V,CP)
% draws the element edges

mu = length(U);
mv = length(V);
nu = length(CP(:,1,1));
nv = length(CP(1,:,1));

grid=100;
eps=10e-10;

% I) edges in u-direction
P = zeros(31,mu-2*p+mv-2*q);
l=1;  % counting index of lines
for j2 = q+1:mv-q
  v = V(j2);
  j = findspan(v,V,nv);
  s = (U(mu)-U(1))/grid;  %incremental step for u
  u = U(1);
  k=1;
  while u <= U(mu)+eps
    i = findspan(u,U,nu);
    P(k,l,1:3) = get_point_surf(p,i,u,U,q,j,v,V,CP);
    k=k+1;
    u=u+s;
  end
  l=l+1;
end

% II) edges in v-direction
for i2 = p+1:mu-p
  u = U(i2);
  i = findspan(u,U,nu);
  s = (V(mv)-V(1))/grid;  %incremental step for v
  v = V(1);
  k=1;
  while v <= V(mv)+eps
    j = findspan(v,V,nv);
    P(k,l,1:3) = get_point_surf(p,i,u,U,q,j,v,V,CP);
    k=k+1;
    v=v+s;
  end
  l=l+1;
end

plot3 (P(:,:,1),P(:,:,2),P(:,:,3),'Color','black');