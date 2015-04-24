function rhs = rhs_grd_sigma_ANS(p,U,q,V,r,W,CP,d,E,nue,ngauss)
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

% initialize global vector
rhs  = zeros(ndof,1);

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
        rhs_e = zeros(ndof_e,1);

        % loops over integration points
        [GPw,GWw] = gauss(ngauss(3));
        for kw = 1:ngauss(3)
          
          % all ANS operations inside zeta loop
          Ba = zeros(6,ndof_e,9);
          %S = zeros(3,3,27);
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
              Bc(:,:,itp) = B_matrix_curv(i,p,u,U,j,q,v,V,k,r,w,W,CP);
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
              Bc(:,:,itp) = B_matrix_curv(i,p,u,U,j,q,v,V,k,r,w,W,CP);
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
              Bc(:,:,itp) = B_matrix_curv(i,p,u,U,j,q,v,V,k,r,w,W,CP);
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
              temp = B_matrix_curv(i,p,u,U,j,q,v,V,k,r,w,W,CP);
              Ba(3,:,iip) = temp(3,:);
              iip = iip+1;
            end
          end
          
          % --------------------------------------
          % Then: Loop over integration points:
          iip = 1;
          for kv = 1:ngauss(2)
            for ku = 1:ngauss(1)
              % NURBS coordinates u,v from gauss coordinates
              u = ( U(i+1)+U(i) + GPu(ku)*(U(i+1)-U(i)) )/2;
              v = ( V(j+1)+V(j) + GPv(kv)*(V(j+1)-V(j)) )/2;
              w = ( W(k+1)+W(k) + GPw(kw)*(W(k+1)-W(k)) )/2;
              dJ = detJ(i,p,u,U,j,q,v,V,k,r,w,W,CP);
              S = get_stress_curv(p,q,r,u,v,w,U,V,W,CP,d,E,nue);
              S_vet = [S(1,1); S(2,2); S(3,3); S(1,2); S(2,3); S(3,1)];
              gwu = GWu(ku); gwv = GWv(kv); gww = GWw(kw);
              rhs_e = Ba(:,:,iip)'*S_vet*gwu*gwv*gww*dJ*map + rhs_e;
              iip = iip+1;
            end
          end
          
        end
        
        % Insert rhs_e into global vector rhs
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
         % insert rhs_e into rhs
        for row = 1:3*(1+p)*(1+q)*(1+r)
            rhs(dof_l(row),1) = rhs_e(row,1) + rhs(dof_l(row),1);
        end
      end
    end
  end
end