function S = get_point_solid(u,v,w,i,j,k,p,q,r,U,V,W,CP)
% returns X,Y,Z coordinates corresponding to u,v,w
% i,j,k are the knot span indices. if unknown, insert 0!

if (i==0); i = findspan(u,U,length(CP(:,1,1,1)));  end
if (j==0); j = findspan(v,V,length(CP(1,:,1,1)));  end
if (k==0); k = findspan(w,W,length(CP(1,1,:,1)));  end

Nu=basisfunc(i,p,u,U);
Nv=basisfunc(j,q,v,V);
Nw=basisfunc(k,r,w,W);

SumNw = 0;
for c = 0:r
  for b = 0:q 
    for a = 0:p
      SumNw = Nu(a+1)*Nv(b+1)*Nw(c+1)*CP(i-p+a,j-q+b,k-r+c,4)+SumNw;
    end
  end
end

S(1:3) = 0;
for c = 0:r
  for b = 0:q 
    for a = 0:p
      S(1) = Nu(a+1)*Nv(b+1)*Nw(c+1)*CP(i-p+a,j-q+b,k-r+c,4)*CP(i-p+a,j-q+b,k-r+c,1)/SumNw+S(1);
      S(2) = Nu(a+1)*Nv(b+1)*Nw(c+1)*CP(i-p+a,j-q+b,k-r+c,4)*CP(i-p+a,j-q+b,k-r+c,2)/SumNw+S(2);
      S(3) = Nu(a+1)*Nv(b+1)*Nw(c+1)*CP(i-p+a,j-q+b,k-r+c,4)*CP(i-p+a,j-q+b,k-r+c,3)/SumNw+S(3); 
    end
  end
end