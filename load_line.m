function fl = load_line(fl_old,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,fload,dir,proj)
% returns the consistent nodal forces to a line load f.
% direction of f can be in x,y,z 
% Parameters:
%     fl_old      existing force vector
%     fl:         updated force vector
%     ub,vb,wb:   load extension (e.g. ub=[0 1], vb=[0 0], wb=[1 1])
%     ngauss:     vector with no. of gauss points for u, v and w
%     fload:      constant line load or handle to load function in [N/m]
%     dir:        load direction: 1-x, 2-y, 3-z
%     proj:       flag if load must be projected onto line
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
  mapu = 1;
else
  i2 = findspan(ub(2),U,nu);
  if (ub(2)~=U(mu));    i2=i2-1;    end
  ugauss = ngauss(1);   uvw = 1;
end
%---------------- V ---------------------
j1 = findspan(vb(1),V,nv);
if (isscalar(vb))
  v = vb;       j2 = j1;
  vgauss = 1;   gwv = 1;
  mapv = 1;
else
  j2 = findspan(vb(2),V,nv);
  if (vb(2)~=V(mv));    j2=j2-1;    end
  vgauss = ngauss(2);   uvw = 2;
end
%---------------- W ---------------------
k1 = findspan(wb(1),W,nw);
if (isscalar(wb))
  w = wb;       k2 = k1;
  wgauss = 1;   gww = 1;
  mapw = 1;
else
  k2 = findspan(wb(2),W,nw);
  if (wb(2)~=W(mw));    k2=k2-1;    end
  wgauss = ngauss(3);   uvw = 3;
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
              if (isnumeric(fload)==0); f=fload(p,q,r,i,j,k,u,v,w,U,V,W,CP); end
              ds = integral_edge(p,q,r,i,j,k,u,v,w,U,V,W,CP,uvw,proj);
              if     (proj==0);
                Fel(:,:,:,dir) = f*map*gw*ds;
              elseif (proj==1);  % projection of ds to the plane perpendicular to dir
                e1 = dir+1;  if(e1>3); e1=e1-3; end;  ds_prj1=ds(:,:,:,e1).^2;
                e2 = dir+2;  if(e2>3); e2=e2-3; end;  ds_prj2=ds(:,:,:,e2).^2;
                ds_prj = sqrt(ds_prj1 + ds_prj2);
                Fel(:,:,:,dir) = f*map*gw*ds_prj;
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