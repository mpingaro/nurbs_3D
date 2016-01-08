function K = stiff_mat(p,U,q,V,r,W,CP,E,nue,ngauss)
% Returns global stiffness matrix for a plate in plane stress element
% Input: p,q,r:    polynomial degrees
%        U,V,W:    knot vectors
%        CP:       control points
%        E,nue:    material parameter
%        ngauss:   vector with no. of gauss points for u, v and w

% Read input
 % Control Points CP
mu = length(U);
mv = length(V);
mw = length(W);
nu = length(CP(:,1,1,1));
nv = length(CP(1,:,1,1));
nw = length(CP(1,1,:,1));
check_input(p,mu,nu,q,mv,nv,r,mw,nw);
 % degrees of freedom ndof
ndof = 3*nu*nv*nw;

 % material matrix D
D = E/((1+nue)*(1-2*nue))*[1-nue nue nue 0 0 0; nue 1-nue nue 0 0 0; nue nue 1-nue 0 0 0;...
    0 0 0 (1-2*nue)/2 0 0; 0 0 0 0 (1-2*nue)/2 0; 0 0 0 0 0 (1-2*nue)/2];

% dof assigns three DoF (x,y,z) to every CP
%  numbering follows CP: CP1->dof1,dof2,dof3 CP2->dof4,dof5,dof6
dof = zeros(nu,nv,nw,3);
k=1;
for cpk = 1:nw
  for cpj = 1:nv
    for cpi = 1:nu
      dof(cpi,cpj,cpk,1)=k;
      dof(cpi,cpj,cpk,2)=k+1;
      dof(cpi,cpj,cpk,3)=k+2;
      k=k+3;
    end
  end
end

% initialize global stiffness matrix
K  = zeros(ndof,ndof);

% loops over elements
for k = (r+1):(mw-r-1)
  for j = (q+1):(mv-q-1)
    for i = (p+1):(mu-p-1)
      % check if element is greater than zero
      if (U(i+1)~=U(i) && V(j+1)~=V(j) && W(k+1)~=W(k))
        ndof_e = 3*(p+1)*(q+1)*(r+1);
        ke = zeros(ndof_e,ndof_e);

        % Element stiffness matrix, loop over Gauss points
        [GPu,GWu] = gauss(ngauss(1));
        [GPv,GWv] = gauss(ngauss(2));
        [GPw,GWw] = gauss(ngauss(3));
        for kw = 1:ngauss(3)
          for kv = 1:ngauss(2)
            for ku = 1:ngauss(1)
              % NURBS coordinates u,v from gauss coordinates
              u = ( U(i+1)+U(i) + GPu(ku)*(U(i+1)-U(i)) )/2;
              v = ( V(j+1)+V(j) + GPv(kv)*(V(j+1)-V(j)) )/2;
              w = ( W(k+1)+W(k) + GPw(kw)*(W(k+1)-W(k)) )/2;
              gwu = GWu(ku);
              gwv = GWv(kv);
              gww = GWw(kw);
              
%               % material matrix D
%               F0 = jacobian(i,p,u,U,j,q,v,V,k,r,w,W,CP0);
%               Ft = jacobian(i,p,u,U,j,q,v,V,k,r,w,W,CP);
%               %F = Ft*inv(F0);
%               F =  Ft*(F0\blkdiag(1,1,1));
%               D = material_elastic_matrix(E, nue, F);
              
              ke = stiff_mat_el(i,p,u,U,j,q,v,V,k,r,w,W,CP,gwu,gwv,gww,D) + ke;
            end
          end
        end
        
        % Insert ke into global stiffness matrix K
         % relation global-local dof
        dof_l = zeros(1,ndof_e);
        l=1;
        for cpk = k-r:k
          for cpj = j-q:j
            for cpi = i-p:i
              dof_l(l)  =dof(cpi,cpj,cpk,1);
              dof_l(l+1)=dof(cpi,cpj,cpk,2);
              dof_l(l+2)=dof(cpi,cpj,cpk,3);
              l=l+3;
            end
          end
        end
         % insert ke into K
        for col = 1:3*(1+p)*(1+q)*(1+r)
          for row = 1:3*(1+p)*(1+q)*(1+r)
            K(dof_l(row),dof_l(col)) = ke(row,col) + K(dof_l(row),dof_l(col));
          end
        end
      end
    end
  end
end