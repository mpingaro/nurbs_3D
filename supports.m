function rb = supports(rb,u,v,w,dir,CP)
% makes fixed supports and adds these to existing ones
% suitable to fix a corner, an edge or a face 
% Parameters:
%     rb:      supports (input and output)
%     u,v,w:   region to be supported (e.g. u=[0 1], v=[0 1], w=[0 0])
%     dir:     direction: 1-x, 2-y, 3-z

r = length(rb)+1;
nu = length(CP(:,1,1,1));
nv = length(CP(1,:,1,1));
nw = length(CP(1,1,:,1));

for k = w(1)*(nw-1)+1:w(2)*(nw-1)+1
  for j = v(1)*(nv-1)+1:v(2)*(nv-1)+1
    for i = u(1)*(nu-1)+1:u(2)*(nu-1)+1
      rb(r) = 3*((k-1)*nu*nv + (j-1)*nu + i-1)+dir;
      r=r+1;
    end
  end
end

% sort rb and delete double entries
rb = sort(rb);
i=1;
while i < length(rb)
  if (rb(i)==rb(i+1));  rb(i+1)=[];  i=i-1;  end
  i=i+1;
end