function fl = load_body(fl_old,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,fload,dir)
% returns the consistent nodal forces to a body load fload.
% Parameters:
%     fl_old      existing force vector
%     fl:         updated force vector
%     ub,vb,wb:   load extension (e.g. ub=[0 0.5], vb=[0 1], wb=[0 1])
%     ngauss:     vector with no. of gauss points for u, v and w
%     fload:      constant body load or handle to load function in [N/m^3]
%     dir:        load direction: 1-x, 2-y, 3-z

mu = length(U);
mv = length(V);
mw = length(W);
nu = length(CP(:,1,1,1));
nv = length(CP(1,:,1,1));
nw = length(CP(1,1,:,1));

i1 = findspan(ub(1),U,nu);
i2 = findspan(ub(2),U,nu);  if (ub(2)~=U(mu));  i2=i2-1;  end
j1 = findspan(vb(1),V,nv);
j2 = findspan(vb(2),V,nv);  if (vb(2)~=V(mv));  j2=j2-1;  end
k1 = findspan(wb(1),W,nw);
k2 = findspan(wb(2),W,nw);  if (wb(2)~=W(mw));  k2=k2-1;  end

if (isnumeric(fload)==1);  f = fload;  end
% ------------------------- loop over elements ----------------------------
F = zeros(nu,nv,nw,3);
for k = k1:k2
  for j = j1:j2
    for i = i1:i2
      % check if element is greater than zero
      if (U(i+1)~=U(i) && V(j+1)~=V(j) && W(k+1)~=W(k))
        % ---------- loop over Gauss points -------------------------------
        [GPu,GWu] = gauss(ngauss(1));
        [GPv,GWv] = gauss(ngauss(2));
        [GPw,GWw] = gauss(ngauss(3));
        for kw = 1:ngauss(3)
          for kv = 1:ngauss(2)
            for ku = 1:ngauss(1)
              u = ( U(i+1)+U(i) + GPu(ku)*(U(i+1)-U(i)) )/2;
              mapu = (U(i+1)-U(i))/2;
              gwu = GWu(ku);
              v = ( V(j+1)+V(j) + GPv(kv)*(V(j+1)-V(j)) )/2;
              mapv = (V(j+1)-V(j))/2;
              gwv = GWv(kv);
              w = ( W(k+1)+W(k) + GPw(kw)*(W(k+1)-W(k)) )/2;
              mapw = (W(k+1)-W(k))/2;
              gww = GWw(kw);
              map = mapu*mapv*mapw;
              gw = gwu*gwv*gww;
              if (isnumeric(fload)==0); f=fload(p,q,r,i,j,k,u,v,w,U,V,W,CP); end
              dV = integral_volume(p,q,r,i,j,k,u,v,w,U,V,W,CP);
              Fel = f*map*gw*dV;
              F(i-p:i,j-q:j,k-r:k,dir) = Fel(:,:,:)+F(i-p:i,j-q:j,k-r:k,dir);
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