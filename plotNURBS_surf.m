function plotNURBS_surf(p,q,U,V,CP)
% draws a NURBS surface grid 50x50 lines

mu = length(U);
mv = length(V);
nu = length(CP(:,1,1));
nv = length(CP(1,:,1));

grid=30;
eps=10e-10;

l=1;  % counting index of lines
r=(V(mv)-V(1))/grid;    %incremental step for v
P = zeros(31,31,3);
v=V(1);
while v <= V(mv)+eps
  j = findspan(v,V,nv);
  s=(U(mu)-U(1))/grid;  %incremental step for u
  u=U(1);
  k=1;
  while u <= U(mu)+eps
    i = findspan(u,U,nu);
    P(k,l,1:3) = get_point_surf(p,i,u,U,q,j,v,V,CP);
    k=k+1;
    u=u+s;
  end
  l=l+1;
  v=v+r;
end

surf(P(:,:,1),P(:,:,2),P(:,:,3),'FaceColor','green','EdgeColor','none');
% camlight left; lighting phong;
% axis equal;
% xlabel('x','FontSize',14);
% ylabel('y','FontSize',14);
% zlabel('z','FontSize',14);