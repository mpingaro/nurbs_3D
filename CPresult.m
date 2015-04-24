function CPd = CPresult(CP,d)
% takes the undeformed control point coordinates and the displacement vector and returns
% the deformed control points

a=1;
CPd = CP;
for k = 1:length(CP(1,1,:,1))
  for j = 1:length(CP(1,:,1,1))
    for i = 1:length(CP(:,1,1,1))
      CPd(i,j,k,1) = CP(i,j,k,1) + d(a);
      CPd(i,j,k,2) = CP(i,j,k,2) + d(a+1);
      CPd(i,j,k,3) = CP(i,j,k,3) + d(a+2);
      a=a+3;
    end
  end
end