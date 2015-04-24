function K = stiff_mat_ANS_xyz(p,U,q,V,r,W,CP,E,nue,ngauss)
% Returns global stiffness matrix for a plate in plane stress element
% Input: p,q,r:    polynomial degrees
%        U,V,W:    knot vectors
%        CP:       control points
%        E,nue:    material parameter
%        ngauss:   vector with no. of gauss points for u, v and w
% J. Kiendl

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
ndof_e = 3*(p+1)*(q+1)*(r+1);
 % material matrix D
D = E/((1+nue)*(1-2*nue))*[1-nue nue nue 0 0 0; nue 1-nue nue 0 0 0; nue nue 1-nue 0 0 0
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

% for ANS: extrapolation matrices
% 1. for xixi and xizeta
M1 = [1/2+sqrt(9/20)  1/2-sqrt(9/20) 0 0 0 0;...
      1/2             1/2            0 0 0 0;...
      1/2-sqrt(9/20)  1/2+sqrt(9/20) 0 0 0 0;...
      0 0 1/2+sqrt(9/20)  1/2-sqrt(9/20) 0 0;...
      0 0 1/2             1/2            0 0;...
      0 0 1/2-sqrt(9/20)  1/2+sqrt(9/20) 0 0;...
      0 0 0 0 1/2+sqrt(9/20)  1/2-sqrt(9/20);...
      0 0 0 0 1/2             1/2           ;...
      0 0 0 0 1/2-sqrt(9/20)  1/2+sqrt(9/20)];
% 2. for etaeta and etazeta
M2 = [1/2+sqrt(9/20) 0 0 1/2-sqrt(9/20) 0 0;...
      0 1/2+sqrt(9/20) 0 0 1/2-sqrt(9/20) 0;...
      0 0 1/2+sqrt(9/20) 0 0 1/2-sqrt(9/20);...
      1/2            0 0 1/2            0 0;...
      0 1/2            0 0 1/2            0;...
      0 0 1/2            0 0 1/2           ;...
      1/2-sqrt(9/20) 0 0 1/2+sqrt(9/20) 0 0;...
      0 1/2-sqrt(9/20) 0 0 1/2+sqrt(9/20) 0;...
      0 0 1/2-sqrt(9/20) 0 0 1/2+sqrt(9/20)];
% 3. for xieta
M3 = 1/4*[14/5+sqrt(36/5) -4/5          -4/5          14/5-sqrt(36/5);...
          1+sqrt(9/5)   1+sqrt(9/5)   1-sqrt(9/5)   1-sqrt(9/5);...
          -4/5          14/5+sqrt(36/5) 14/5-sqrt(36/5) -4/5;...
          1+sqrt(9/5)   1-sqrt(9/5)   1+sqrt(9/5)   1-sqrt(9/5);...
          1           1           1           1;...
          1-sqrt(9/5)   1+sqrt(9/5)   1-sqrt(9/5)   1+sqrt(9/5);...
          -4/5          14/5-sqrt(36/5) 14/5+sqrt(36/5) -4/5;...
          1-sqrt(9/5)   1-sqrt(9/5)   1+sqrt(9/5)   1+sqrt(9/5);...
          14/5-sqrt(36/5) -4/5          -4/5          14/5+sqrt(36/5)];

% loops over elements
for k = (r+1):(mw-r-1)
  for j = (q+1):(mv-q-1)
    for i = (p+1):(mu-p-1)
      % check if element is greater than zero
      if (U(i+1)~=U(i) && V(j+1)~=V(j) && W(k+1)~=W(k))
        map = (U(i+1)-U(i))*(V(j+1)-V(j))*(W(k+1)-W(k))/8;
        ke = zeros(ndof_e,ndof_e);

        % loops over integration points
        [GPw,GWw] = gauss(ngauss(3));
        for kw = 1:ngauss(3)
          
          % all ANS operations inside zeta loop
          Ba = zeros(6,ndof_e,9);
          % --------------------------------------
          % 1. ANS for eps_xixi and eps_xizeta
          Bc = zeros(6,ndof_e,6);
          % 1a) loop over tps, compute Bc and save it
          itp = 1;
          for kv = 1:ngauss(2)
            [GPv,GWv] = gauss(ngauss(2));
            for ku = 1:ngauss(1)-1
              [GPu,GWu] = gauss(ngauss(1)-1);
              % NURBS coordinates u,v from gauss coordinates
              u = ( U(i+1)+U(i) + GPu(ku)*(U(i+1)-U(i)) )/2;
              v = ( V(j+1)+V(j) + GPv(kv)*(V(j+1)-V(j)) )/2;
              w = ( W(k+1)+W(k) + GPw(kw)*(W(k+1)-W(k)) )/2;
              Bc(:,:,itp) = B_matrix(i,p,u,U,j,q,v,V,k,r,w,W,CP);
              itp = itp+1;
            end
          end
          % 1b) extrapolating to integration points, only eps_xixi and _xizeta, i.e. B(1,:) and B(5,:)
          Ba(1,:,:) = (M1*(squeeze(Bc(1,:,:)))')';
          Ba(6,:,:) = (M1*(squeeze(Bc(6,:,:)))')';
          % --------------------------------------
          % 2. ANS for eps_etaeta and eps_etazeta
          Bc = zeros(6,ndof_e,6);
          % 2a) loop over tps, compute Bc and save it
          itp = 1;
          for kv = 1:ngauss(2)-1
            [GPv,GWv] = gauss(ngauss(2)-1);
            for ku = 1:ngauss(1)
              [GPu,GWu] = gauss(ngauss(1));
              % NURBS coordinates u,v from gauss coordinates
              u = ( U(i+1)+U(i) + GPu(ku)*(U(i+1)-U(i)) )/2;
              v = ( V(j+1)+V(j) + GPv(kv)*(V(j+1)-V(j)) )/2;
              w = ( W(k+1)+W(k) + GPw(kw)*(W(k+1)-W(k)) )/2;
              Bc(:,:,itp) = B_matrix(i,p,u,U,j,q,v,V,k,r,w,W,CP);
              itp = itp+1;
            end
          end
          % 2b) extrapolating to integration points, only eps_xixi and _xizeta, i.e. B(2,:) and B(4,:)
          Ba(2,:,:) = (M2*(squeeze(Bc(2,:,:)))')';
          Ba(5,:,:) = (M2*(squeeze(Bc(5,:,:)))')';
          % --------------------------------------
          % 3. ANS for eps_xieta
          Bc = zeros(6,ndof_e,4);
          % 3a) loop over tps, compute Bc and save it
          itp = 1;
          for kv = 1:ngauss(2)-1
            [GPv,GWv] = gauss(ngauss(2)-1);
            for ku = 1:ngauss(1)-1
              [GPu,GWu] = gauss(ngauss(1)-1);
              % NURBS coordinates u,v from gauss coordinates
              u = ( U(i+1)+U(i) + GPu(ku)*(U(i+1)-U(i)) )/2;
              v = ( V(j+1)+V(j) + GPv(kv)*(V(j+1)-V(j)) )/2;
              w = ( W(k+1)+W(k) + GPw(kw)*(W(k+1)-W(k)) )/2;
              Bc(:,:,itp) = B_matrix(i,p,u,U,j,q,v,V,k,r,w,W,CP);
              itp = itp+1;
            end
          end
          % 2b) extrapolating to integration points, only eps_xixi and _xizeta, i.e. B(2,:) and B(4,:)
          Ba(4,:,:) = (M3*(squeeze(Bc(4,:,:)))')';
          % --------------------------------------
          % 4. Normal strain for eps_zetazeta
          % 3a) loop over ips, compute Bc and save it
          iip = 1;
          for kv = 1:ngauss(2)
            [GPv,GWv] = gauss(ngauss(2));
            for ku = 1:ngauss(1)
              [GPu,GWu] = gauss(ngauss(1));
              % NURBS coordinates u,v from gauss coordinates
              u = ( U(i+1)+U(i) + GPu(ku)*(U(i+1)-U(i)) )/2;
              v = ( V(j+1)+V(j) + GPv(kv)*(V(j+1)-V(j)) )/2;
              w = ( W(k+1)+W(k) + GPw(kw)*(W(k+1)-W(k)) )/2;
              temp = B_matrix(i,p,u,U,j,q,v,V,k,r,w,W,CP);
              Ba(3,:,iip) = temp(3,:);
%               Ba(1,:,iip) = temp(1,:);
%               Ba(2,:,iip) = temp(2,:);
%               Ba(4,:,iip) = temp(4,:);
%               Ba(5,:,iip) = temp(5,:);
%               Ba(6,:,iip) = temp(6,:);
              iip = iip+1;
            end
          end
          % --------------------------------------
          % Then: Loop over integration points: call B-matrix and perform
          % numerical integration of ke
          iip = 1;
          for kv = 1:ngauss(2)
            for ku = 1:ngauss(1)
              % NURBS coordinates u,v from gauss coordinates
              u = ( U(i+1)+U(i) + GPu(ku)*(U(i+1)-U(i)) )/2;
              v = ( V(j+1)+V(j) + GPv(kv)*(V(j+1)-V(j)) )/2;
              w = ( W(k+1)+W(k) + GPw(kw)*(W(k+1)-W(k)) )/2;
              dJ = detJ(i,p,u,U,j,q,v,V,k,r,w,W,CP);
              gwu = GWu(ku); gwv = GWv(kv); gww = GWw(kw);
              ke = Ba(:,:,iip)'*D*Ba(:,:,iip)*gwu*gwv*gww*dJ*map + ke;
              iip = iip+1;
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