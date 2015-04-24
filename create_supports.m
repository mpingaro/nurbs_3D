function [xs,ys,zs] = create_supports(CP,rb)
% draws triangles for supports

nu = length(CP(:,1,1,1));
nv = length(CP(1,:,1,1));

% scaling factor
up = max(max(max(max(CP))));
lo = min(min(min(min(CP))));
fac = (up-lo)/5;

xs = zeros(length(rb),4);
ys = zeros(length(rb),4);
zs = zeros(length(rb),4);
for l = 1:length(rb)
  % get the corresponding Control Point number p and indices CP(i,j)
  h = (rb(l)+2)/3;
  p = floor(h);
  d = round(3*(h-p))+1;
  k = ceil(p/(nu*nv));
  j = ceil((p-(k-1)*nu*nv)/nu);
  i = p-(k-1)*nu*nv-(j-1)*nu;
  if (d==1)              %(x-support)
    xs(l,1)=CP(i,j,k,1);
    xs(l,2)=CP(i,j,k,1)-0.1732*fac;
    xs(l,3)=CP(i,j,k,1)-0.1732*fac;
    xs(l,4)=xs(l,1);
    ys(l,1)=CP(i,j,k,2);
    ys(l,2)=CP(i,j,k,2)+0.1*fac;
    ys(l,3)=CP(i,j,k,2)-0.1*fac;
    ys(l,4)=ys(l,1);
    zs(l,1:4)=CP(i,j,k,3);
  elseif (d==2)          %(y-support)
    xs(l,1)=CP(i,j,k,1);
    xs(l,2)=CP(i,j,k,1)-0.1*fac;
    xs(l,3)=CP(i,j,k,1)+0.1*fac;
    xs(l,4)=xs(l,1);
    ys(l,1)=CP(i,j,k,2);
    ys(l,2)=CP(i,j,k,2)-0.1732*fac;
    ys(l,3)=CP(i,j,k,2)-0.1732*fac;
    ys(l,4)=ys(l,1);
    zs(l,1:4)=CP(i,j,k,3);
  elseif (d==3)          %(y-support)
    xs(l,1)=CP(i,j,k,1);
    xs(l,2)=CP(i,j,k,1)-0.1*fac;
    xs(l,3)=CP(i,j,k,1)+0.1*fac;
    xs(l,4)=xs(l,1);
    ys(l,1:4)=CP(i,j,k,2);
    zs(l,1)=CP(i,j,k,3);
    zs(l,2)=CP(i,j,k,3)-0.1732*fac;
    zs(l,3)=CP(i,j,k,3)-0.1732*fac;
    zs(l,4)=zs(l,1);
  end
end