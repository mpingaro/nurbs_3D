function [xf,yf,zf] = create_arrows(CP,f)

d=1;
l=1;
for k = 1:length(CP(1,1,:,1))
  for j = 1:length(CP(1,:,1,1))
    for i = 1:length(CP(:,1,1,1))
      if (f(d)~=0)
        xf(l,1) = CP(i,j,k,1);
        xf(l,2) = CP(i,j,k,1)-f(d)/max(abs(f));
        yf(l,1) = CP(i,j,k,2);
        yf(l,2) = CP(i,j,k,2);
        zf(l,1) = CP(i,j,k,3);
        zf(l,2) = CP(i,j,k,3);
        l=l+1;
      end
      if (f(d+1)~=0)
        xf(l,1) = CP(i,j,k,1);
        xf(l,2) = CP(i,j,k,1);
        yf(l,1) = CP(i,j,k,2);
        yf(l,2) = CP(i,j,k,2)-f(d+1)/max(abs(f));
        zf(l,1) = CP(i,j,k,3);
        zf(l,2) = CP(i,j,k,3);
        l=l+1;
      end
      if (f(d+2)~=0)
        xf(l,1) = CP(i,j,k,1);
        xf(l,2) = CP(i,j,k,1);
        yf(l,1) = CP(i,j,k,2);
        yf(l,2) = CP(i,j,k,2);
        zf(l,1) = CP(i,j,k,3);
        zf(l,2) = CP(i,j,k,3)-f(d+2)/max(abs(f));
        l=l+1;
      end
      d=d+3;
    end
  end
end