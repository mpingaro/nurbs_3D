function CP = CP_matrix(p,q,r,U,V,W,CPm)
% writes control points from a list (x*4 matrix) into a nu*nv*nw matrix 
% (nu*nv*nw*4 3D matrix). nu,nv,nw are determined from U,V,W and p,q,r

mu=length(U); nu=mu-p-1;
mv=length(V); nv=mv-q-1;
mw=length(W); nw=mw-r-1;
CP = zeros(nu,nv,nw,4);
l=1;
% for k=1:nw
%   for j=1:nv
%     for i=1:nu
%       CP(i,j,k,:)=CPm(l,:);
%       l=l+1;
%     end
%   end
% end
for i=1:nu
  for j=1:nv
    for k=1:nw
      CP(i,j,k,:)=CPm(l,:);
      l=l+1;
    end
  end
end

