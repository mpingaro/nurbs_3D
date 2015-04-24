function [CPr,Ur,Vr,Wr] = knot_refine_solid(p,q,r,U,V,W,CP,Ru,Rv,Rw)
%
% adapted from algorithm in Piegl, Les. "The NURBS Book". Springer-Verlag: 
% Berlin 1995; p. 167.

nu=length(CP(:,1,1,1));
nv=length(CP(1,:,1,1));
nw=length(CP(1,1,:,1));
%==================== Projective control points Pw ========================
for k = 1:nw
  for j = 1:nv
    for i = 1:nu
      Pw(i,j,k,1:3) = CP(i,j,k,1:3)*CP(i,j,k,4);
      Pw(i,j,k,4)   = CP(i,j,k,4);
    end
  end
end
%============================ Refine in u =================================
if isempty(Ru)
  Ur=U;
  Qw=Pw;
else
  rr = length(Ru); 
  mu = nu+p+1;
  a = findspan(Ru(1),U,nu);
  b = findspan(Ru(rr),U,nu)+1;

  for k = 1:nw
    for j = 1:nv
      for  i = 1:a-p;    Qw(i,j,k,:)    = Pw(i,j,k,:);   end
      for  i = b-1:nu;   Qw(i+rr,j,k,:) = Pw(i,j,k,:);   end
    end
  end
  for i = 1:a;      Ur(i)    = U(i);  end
  for i = b+p:mu;   Ur(i+rr) = U(i);  end
    
  i = b+p-1;   ir = i+rr;
  for  iru = rr:-1:1
    while (Ru(iru)<=U(i) && i>a)
      for k = 1:nw
        for j = 1:nv
          Qw(ir-p-1,j,k,:) = Pw(i-p-1,j,k,:);
        end
      end
      Ur(ir) = U(i);
      ir = ir-1;   i = i-1;
    end
    
    for k = 1:nw
      for j = 1:nv
        Qw(ir-p-1,j,k,:) = Qw(ir-p,j,k,:);
      end
    end
    for l = 1:p
      ind = ir-p+l;
      alfa = (Ru(iru)-Ur(ir+l)) / (U(i-p+l)-Ur(ir+l));
      for k = 1:nw
        for j = 1:nv
          Qw(ind-1,j,k,:) = alfa*Qw(ind-1,j,k,:) + (1-alfa)*Qw(ind,j,k,:);
        end
      end
    end
    Ur(ir) = Ru(iru);
    ir = ir-1;
  end
  Pw = Qw;
  nu = nu+rr;
end
%============================ Refine in v =================================
if isempty(Rv)
  Vr=V;
else
  rr=length(Rv); 
  mv = nv+q+1;
  a = findspan(Rv(1),V,nv);
  b = findspan(Rv(rr),V,nv)+1;
  
  for k = 1:nw
    for i = 1:nu
      for j = 1:a-q;    Qw(i,j,k,:)  = Pw(i,j,k,:);    end
      for j = b-1:nv;   Qw(i,j+rr,k,:) = Pw(i,j,k,:);  end
    end
  end 
  for j = 1:a;      Vr(j)    = V(j);  end
  for j = b+q:mv;   Vr(j+rr) = V(j);  end

  j = b+q-1;   jr = j+rr;
  for jrv = rr:-1:1
    while ((Rv(jrv)<=V(j)) && (j>a))
      for k = 1:nw
        for i = 1:nu
          Qw(i,jr-q-1,k,:) = Pw(i,j-q-1,k,:);
        end
      end
      Vr(jr) = V(j);
      jr = jr-1;   j = j-1;
    end

    for k = 1:nw
      for i = 1:nu
        Qw(i,jr-q-1,k,:) = Qw(i,jr-q,k,:);
      end
    end
    for l = 1:q
      ind = jr-q+l;
      alfa = (Rv(jrv)-Vr(jr+l)) / (V(j-q+l)-Vr(jr+l));
      for k = 1:nw
        for i = 1:nu
          Qw(i,ind-1,k,:) = alfa*Qw(i,ind-1,k,:) + (1-alfa)*Qw(i,ind,k,:);
        end
      end
    end
    Vr(jr) = Rv(jrv);
    jr = jr-1;
  end
  Pw = Qw;
  nv = nv+rr;
end
%============================ Refine in w =================================
if isempty(Rw)
  Wr=W;
else
  rr=length(Rw); 
  mw = nw+r+1;
  a = findspan(Rw(1),W,nw);
  b = findspan(Rw(rr),W,nw)+1;
  
  for j = 1:nv
    for i = 1:nu
      for k = 1:a-r;    Qw(i,j,k,:)  = Pw(i,j,k,:);    end
      for k = b-1:nw;   Qw(i,j,k+rr,:) = Pw(i,j,k,:);  end
    end
  end 
  for k = 1:a;      Wr(k)    = W(k);  end
  for k = b+r:mw;   Wr(k+rr) = W(k);  end

  k = b+r-1;   kr = k+rr;
  for krw = rr:-1:1
    while ((Rw(krw)<=W(k)) && (k>a))
      for j = 1:nv
        for i = 1:nu
          Qw(i,j,kr-r-1,:) = Pw(i,j,k-r-1,:);
        end
      end
      Wr(kr) = W(k);
      kr = kr-1;   k = k-1;
    end

    for j = 1:nv
      for i = 1:nu
        Qw(i,j,kr-r-1,:) = Qw(i,j,kr-r,:);
      end
    end
    for l = 1:r
      ind = kr-r+l;
      alfa = (Rw(krw)-Wr(kr+l)) / (W(k-r+l)-Wr(kr+l));
      for j = 1:nv
        for i = 1:nu
          Qw(i,j,ind-1,:) = alfa*Qw(i,j,ind-1,:) + (1-alfa)*Qw(i,j,ind,:);
        end
      end
    end
    Wr(kr) = Rw(krw);
    kr = kr-1;
  end
end
%==================== Retransform from Qw to CPr ==========================
for k = 1:length(Qw(1,1,:,1))
  for j = 1:length(Qw(1,:,1,1))
    for i = 1:length(Qw(:,1,1,1))
      CPr(i,j,k,1:3) = Qw(i,j,k,1:3)/Qw(i,j,k,4);
      CPr(i,j,k,4)   = Qw(i,j,k,4);
    end
  end
end