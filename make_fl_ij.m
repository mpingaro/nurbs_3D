function FL = make_fl_ij(fl,nu,nv)

k=1;
for j = 1:nv
  for i = 1:nu
    FL(i,j,1)=fl(k);
    FL(i,j,2)=fl(k+1);
    k=k+2;
  end
end