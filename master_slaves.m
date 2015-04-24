function ms = master_slaves(closed,CP)
% writes the master-slave-matrix ms
% 1st column = slave    2nd column = master
% closed: vector with flags for u,v,w;  0-not closed  1-closed      

nu = length(CP(:,1,1,1));
nv = length(CP(1,:,1,1));
nw = length(CP(1,1,:,1));

l=0;
if (closed(1) == 1)
  for k = 1:nw
    for j = 1:nv
      for d = 1:3
        ms(l+d,2) = 3*((k-1)*nu*nv + (j-1)*nu) + d;
        ms(l+d,1) = 3*((k-1)*nu*nv + (j-1)*nu + nu-1) + d;
      end
      l=l+3;
    end
  end
end
if (closed(2) == 1)
  for k = 1:nw
    for i = 1:nu
      for d = 1:3
        ms(l+d,2) = 3*((k-1)*nu*nv + i-1) + d;
        ms(l+d,1) = 3*((k-1)*nu*nv + (nv-1)*nu + i-1) + d;
      end
      l=l+3;
    end
  end
end
if (closed(3) == 1)
  for j = 1:nv
    for i = 1:nu
      for d = 1:3
        ms(l+d,2) = 3*((j-1)*nu + i-1) + d;
        ms(l+d,1) = 3*((nw-1)*nu*nv + (j-1)*nu + i-1) + d;
      end
      l=l+3;
    end
  end
end
ms = sortrows(ms,1);
i=1;
while i < length(ms(:,1))
  if (ms(i,1)==ms(i+1,1));  ms(i+1,:)=[];  i=i-1;  end
  i=i+1;
end