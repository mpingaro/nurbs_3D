function fl = load_point(fl_old,u,v,p,q,U,V,CP,f,dir)
% returns the consistent nodal forces to a point load f.
% Parameters:
%     fl_old      existing force vector
%     fl:         updated force vector
%     u,v:        load position
%     f:          point load
%     dir:        direction of f  1=x, 2=y, 3=z

nu = length(CP(:,1,1));
nv = length(CP(1,:,1));
i = findspan(u,U,nu);
j = findspan(v,V,nv);
Nu=basisfunc(i,p,u,U);
Nv=basisfunc(j,q,v,V);

SumNw = 0;
for c = 0:q
  for b = 0:p
    SumNw = Nu(b+1)*Nv(c+1)*CP(i-p+b,j-q+c,4)+SumNw;
  end
end

FL = zeros(nu,nv,2);
for c = 0:q 
  for b = 0:p
    FL(i-p+b,j-q+c,dir) = Nu(b+1)*Nv(c+1)*CP(i-p+b,j-q+c,4)*f/SumNw;
  end
end

fl = make_fl_dof(FL);
if isvector(fl_old);  fl = fl + fl_old;   end