function fl = load_area(fl_old,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,fload,dir,proj)
% returns the consistent nodal forces to an area load f.
% direction of f can be in x,y,z,u,v,w or perpendicular to the face 
% Parameters:
%     fl_old      existing force vector
%     fl:         updated force vector
%     ub,vb,wb:   load extension (e.g. ub=[0 1], vb=[0 1], wb=1)
%     ngauss:     vector with no. of gauss points for u, v and w
%     fload:      constant area load or handle to load function in [N/m^2]
%     dir:        load direction: 1-x, 2-y, 3-z
%                                 4-u  5-v  6-w
%                                 7-perpendicular to face
%     proj:       flag if load must be projected onto face (only for dir=1-3)
%                    0 no projection e.g. self weight
%                    1 load must be projected e.g. snow, wind

mu = length(U);
mv = length(V);
mw = length(W);
nu = length(CP(:,1,1,1));
nv = length(CP(1,:,1,1));
nw = length(CP(1,1,:,1));

if (isnumeric(fload)==1);  f = fload;  end
% ================ Set up parameters for integration ======================
%---------------- U ---------------------
i1 = findspan(ub(1),U,nu);
if (isscalar(ub))
  u = ub;       i2 = i1;
  ugauss = 1;   gwu = 1;
  mapu = 1;     uvw = 1;
else
  i2 = findspan(ub(2),U,nu);
  if (ub(2)~=U(mu));    i2=i2-1;    end
  ugauss = ngauss(1);
end
%---------------- V ---------------------
j1 = findspan(vb(1),V,nv);
if (isscalar(vb))
  v = vb;       j2 = j1;
  vgauss = 1;   gwv = 1;
  mapv = 1;     uvw = 2;
else
  j2 = findspan(vb(2),V,nv);
  if (vb(2)~=V(mv));    j2=j2-1;    end
  vgauss = ngauss(2);
end
%---------------- W ---------------------
k1 = findspan(wb(1),W,nw);
if (isscalar(wb))
  w = wb;       k2 = k1;
  wgauss = 1;   gww = 1;
  mapw = 1;     uvw = 3;
else
  k2 = findspan(wb(2),W,nw);
  if (wb(2)~=W(mw));    k2=k2-1;    end
  wgauss = ngauss(3);
end
% =========================================================================

% ============================ Integration ================================
F = zeros(nu,nv,nw,3);   Fel = zeros(p+1,q+1,r+1,3);
% ------------------------- loop over elements ----------------------------
for k = k1:k2
  for j = j1:j2
    for i = i1:i2
      % check if element is greater than zero
      if (U(i+1)~=U(i) && V(j+1)~=V(j) && W(k+1)~=W(k))
        % ---------- loop over Gauss points -------------------------------
        [GPu,GWu] = gauss(ngauss(1));
        [GPv,GWv] = gauss(ngauss(2));
        [GPw,GWw] = gauss(ngauss(3));
        for kw = 1:wgauss
          for kv = 1:vgauss
            for ku = 1:ugauss
              if (isscalar(ub)==0)
                u = ( U(i+1)+U(i) + GPu(ku)*(U(i+1)-U(i)) )/2;
                mapu = (U(i+1)-U(i))/2;
                gwu = GWu(ku);
              end
              if (isscalar(vb)==0)
                v = ( V(j+1)+V(j) + GPv(kv)*(V(j+1)-V(j)) )/2;
                mapv = (V(j+1)-V(j))/2;
                gwv = GWv(kv);
              end
              if (isscalar(wb)==0)
                w = ( W(k+1)+W(k) + GPw(kw)*(W(k+1)-W(k)) )/2;
                mapw = (W(k+1)-W(k))/2;
                gww = GWw(kw);
              end
              map = mapu*mapv*mapw;
              gw = gwu*gwv*gww;
              if (isnumeric(fload)==0); f=fload(p,q,r,i,j,k,u,v,w,U,V,W,CP);  end
              if (dir==1 || dir==2 || dir==3);
                dA = integral_area(p,q,r,i,j,k,u,v,w,U,V,W,CP,uvw,proj);
                if (proj==0);      Fel(:,:,:,dir) = f*map*gw*dA;
                elseif (proj==1);  Fel(:,:,:,dir) = f*map*gw*dA(:,:,:,dir);
                end
              elseif (dir==4 || dir==5 || dir==6);
                % direction of f is obtained from ds
                dA = integral_area(p,q,r,i,j,k,u,v,w,U,V,W,CP,uvw,0);
                ds = integral_edge_single(p,q,r,i,j,k,u,v,w,U,V,W,CP,dir-3,1);
                norm_ds = sqrt(ds'*ds);
                for e=1:3;   Fel(:,:,:,e) = f*map*gw*dA*ds(e)/norm_ds;  end
              elseif (dir==7);
                % f*projections of dA to the three basic planes
                dA = integral_area(p,q,r,i,j,k,u,v,w,U,V,W,CP,uvw,1);
                Fel(:,:,:,:) = f*map*gw*dA(:,:,:,:);
              end
              F(i-p:i,j-q:j,k-r:k,:) = Fel(:,:,:,:)+F(i-p:i,j-q:j,k-r:k,:);
            end
          end
        end % ------------------------- end loop Gauss points -------------
      end
    end
  end
end % ------------------------- end loop elements -------------------------
% =========================================================================

fl = make_fl_dof(F);
if isvector(fl_old);  fl = fl + fl_old;   end