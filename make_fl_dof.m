function fl = make_fl_dof(FL)

nu = length(FL(:,1,1,1));
nv = length(FL(1,:,1,1));
nw = length(FL(1,1,:,1));
fl = zeros(3*nu*nv*nw,1);
l=1;
for k = 1:nw
  for j = 1:nv
    for i = 1:nu
      fl(l)   = FL(i,j,k,1);
      fl(l+1) = FL(i,j,k,2);
      fl(l+2) = FL(i,j,k,3);
      l=l+3;
    end
  end
end