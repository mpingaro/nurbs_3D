function [pr,qr,rr,Ur,Vr,Wr,CPr] = degree_elevate_solid(p,q,r,U,V,W,CP,tp,tq,tr)
% degree elevate p,q and r by tp,tq and tr

% adapted from algorithm in Piegl, Les. "The NURBS Book". Springer-Verlag:
% Berlin 1995; p. 206.

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
%========================== DEGREE ELEVATE p ==============================
if (tp>0)
  mu = nu+p+1;
  pr = p+tp;
  % preallocate Ur
  ku = 1;
  for i=1:mu-1;  if(U(i)~=U(i+1));  ku = ku+1;   end;   end 
  Ur = zeros(1,mu+tp*ku);
  % compute Bezier degree elevation coefficients
  bezalfs(1,1) = 1.0;
  bezalfs(pr+1,p+1) = 1.0;

  for i = 1:fix(pr/2)   % bezalfs are symmetric to row pr/2
    inv = 1/nchoosek(pr,i);
    mpi = min(p,i);
    for j = max(0,i-tp):mpi   % alfas for degree elevation
      bezalfs(i+1,j+1) = inv*nchoosek(p,j)*nchoosek(tp,i-j);
      bezalfs(pr-i+1,p-j+1) = bezalfs(i+1,j+1);    % symmetry !
    end
  end

  mh = pr;
  kind = pr+1;
  in=-1;
  a=p+1;
  b=p+2;
  cind=1;
  ua = U(a);
  Qw(1,:,:,:) = Pw(1,:,:,:);

  % left end of Ur
  for i = 1:pr+1
    Ur(i) = ua;
  end

  % initialize first Bezier seg
  for i = 1:p+1
    bpts(i,:,:,:) = Pw(i,:,:,:);
  end

  % big loop through knot vector
  while (b<mu)
    i = b;
    while (b<mu && U(b)==U(b+1))
      b = b+1;      % make b the rightmost ocurrence of ub
    end
    mult = b-i+1;   % multiplicity of ub
    mh = mh+mult+tp;
    ub = U(b);
    oldin = in;       % in from last segment
    in = p-mult;      % ub to be inserted in times
  
    if (oldin>0);   lbz = fix((oldin+2)/2);
    else            lbz = 1;
    end
  
    if (in>0);      rbz = pr-fix((in+1)/2);
    else            rbz = pr;
    end
  
    % insert knots to get Bezier segment
    if (in>0)
      numer = ub - ua;
      alfs = zeros(1,p-mult);
      for k = p:-1:mult+1
       alfs(k-mult) = numer/(U(a+k)-ua);   % alfa for knot insertion
      end
      for j = 1:in   % in times knot insertion
        save = in-j+1;
        s = mult+j;
        for k = p+1:-1:s+1    % new CP due to knot insertion
          bpts(k,:,:) = alfs(k-s)*bpts(k,:,:)+(1-alfs(k-s))*bpts(k-1,:,:);
        end
        Nextbpts(save,:,:,:) = bpts(p+1,:,:,:);
      end
    end
    % end of inserting knots
  
    % degree elevate Bezier
    for i = lbz:pr
      % only points lbz,..,pr are used below
      ebpts(i+1,1:nv,1:nw,1:4) = 0;
      mpi = min(p,i);
      for j = max(0,i-tp):mpi   % new CP due to degree elevation
        ebpts(i+1,:,:,:) = ebpts(i+1,:,:,:) + bezalfs(i+1,j+1)*bpts(j+1,:,:,:);
      end
    end
    % end of degree elevating Bezier
  
    % knot removal ua oldin times
    if (oldin>1)
      first = kind-2;
      last = kind;
      den = ub-ua;
      bet = (ub-Ur(kind))/den;
      for tr = 1:(oldin-1)
        i = first;
        j = last;
        kj = j-kind+1;
        while ((j-i)>tr)
          % loop and compute the new CP for one removal step
          if (i<cind)
            alf = (ub-Ur(i+1))/(ua-Ur(i+1));
            Qw(i+1,:,:,:) = (alf*Qw(i+1,:,:,:)+(1.0-alf)*Qw(i,:,:,:));
          end
          if (j>=lbz)
            if ((j-tr)<=(kind-pr+oldin))
              gam = (ub-Ur(j-tr+1))/den;
              ebpts(kj+1,:,:,:) = gam*ebpts(kj+1,:,:,:) + (1.0-gam)*ebpts(kj+2,:,:,:);
            else
              ebpts(kj+1,:,:,:) = bet*ebpts(kj+1,:,:,:) + (1.0-bet)*ebpts(kj+2,:,:,:);
            end
          end
          i = i+1;
          j = j-1;
          kj = kj-1;
        end
        first = first-1;
        last = last+1;
      end
    end   % end of removing knot, u=U(a)
     
    if (a~=p+1)   % load the knot ua
      for i = 0:(pr-oldin-1)
         Ur(kind+1) = ua;
         kind = kind+1;
      end
    end
  
    for j = lbz:rbz   % load CPs into Qw
      Qw(cind+1,:,:,:) =  ebpts(j+1,:,:,:);
      cind = cind +1;
    end
  
    if (b<mu)   % setup for the next pass through loop
      for j = 0:in-1
        bpts(j+1,:,:,:) = Nextbpts(j+1,:,:,:);
      end
      for j = in:p
        bpts(j+1,:,:,:) = Pw(b-p+j,:,:,:);
      end
      a = b;
      b = b+1;
      ua = ub;
    else   % end knots
      for i = 0:pr
        Ur(kind+i+1) = ub;
      end
    end
  end
  % end of big loop through knot vector
  
elseif (tp==0)
  Ur=U;  pr=p;
  Qw=Pw;
end
%========================== DEGREE ELEVATE q ==============================
if (tq>0)
  clear ('Pw');
  Pw(:,:,:,:)=Qw(:,:,:,:);
  nu = length(Pw(:,1,1,1));
  clear ('Qw');  
  clear ('bpts');
  clear ('Nextbpts');
  clear ('ebpts');
  clear ('bezalfs');
  
  mv = nv+q+1;
  qr = q+tq;
  % preallocate Vr
  kv = 1;
  for i=1:mv-1;  if(V(i)~=V(i+1));  kv = kv+1;   end;   end 
  Vr = zeros(1,mv+tq*kv);

  % compute Bezier degree elevation coefficients
  bezalfs(1,1) = 1.0;
  bezalfs(qr+1,q+1) = 1.0;

  for i = 1:fix(qr/2)   % bezalfs are symmetric to row qr/2
    inv = 1/nchoosek(qr,i);
    mpi = min(q,i);
    for j = max(0,i-tq):mpi   % alfas for degree elevation
      bezalfs(i+1,j+1) = inv*nchoosek(q,j)*nchoosek(tq,i-j);
      bezalfs(qr-i+1,q-j+1) = bezalfs(i+1,j+1);    % symmetry !
    end
  end

  mh = qr;
  kind = qr+1;
  in=-1;
  a=q+1;
  b=q+2;
  cind=1;
  ua = V(a);
  Qw(:,1,:,:) = Pw(:,1,:,:);

  % left end of Vr
  for i = 1:qr+1
    Vr(i) = ua;
  end

  % initialize first Bezier seg
  for i = 1:q+1
    bpts(:,i,:,:) = Pw(:,i,:,:);
  end

  % big loop through knot vector
  while (b<mv)
    i = b;
    while (b<mv && V(b)==V(b+1))
      b = b+1;      % make b the rightmost ocurrence of ub
    end
    mult = b-i+1;   % multiplicity of ub
    mh = mh+mult+tq;
    ub = V(b);
    oldin = in;       % in from last segment
    in = q-mult;     % ub to be inserted in times
  
    if (oldin>0);   lbz = fix((oldin+2)/2);
    else           lbz = 1;
    end
  
    if (in>0);      rbz = qr-fix((in+1)/2);
    else           rbz = qr;
    end
  
    % insert knots to get Bezier segment
    if (in>0)
      numer = ub - ua;
      alfs = zeros(1,q-mult);
      for k = q:-1:mult+1
       alfs(k-mult) = numer/(V(a+k)-ua);   % alfa for knot insertion
      end
      for j = 1:in   % in times knot insertion
        save = in-j+1;
        s = mult+j;
        for k = q+1:-1:s+1    % new CP due to knot insertion
          bpts(:,k,:,:) = alfs(k-s)*bpts(:,k,:,:)+(1-alfs(k-s))*bpts(:,k-1,:,:);
        end
        Nextbpts(:,save,:,:) = bpts(:,q+1,:,:);
      end
    end
    % end of inserting knots
  
    % degree elevate Bezier
    for i = lbz:qr
      % only points lbz,..,qr are used below
      ebpts(1:nu,i+1,1:nw,1:4) = 0;
      mpi = min(q,i);
      for j = max(0,i-tq):mpi   % new CP due to degree elevation
        ebpts(:,i+1,:,:) = ebpts(:,i+1,:,:) + bezalfs(i+1,j+1)*bpts(:,j+1,:,:);
      end
    end
    % end of degree elevating Bezier
  
    % knot removal ua oldin times
    if (oldin>1)
      first = kind-2;
      last = kind;
      den = ub-ua;
      bet = (ub-Vr(kind))/den;
      for tr = 1:(oldin-1)
        i = first;
        j = last;
        kj = j-kind+1;
        while ((j-i)>tr)
          % loop and compute the new CP for one removal step
          if (i<cind)
            alf = (ub-Vr(i+1))/(ua-Vr(i+1));
            Qw(:,i+1,:,:) = (alf*Qw(:,i+1,:,:)+(1.0-alf)*Qw(:,i,:,:));
          end
          if (j>=lbz)
            if ((j-tr)<=(kind-qr+oldin))
              gam = (ub-Vr(j-tr+1))/den;
              ebpts(:,kj+1,:,:) = gam*ebpts(:,kj+1,:,:) + (1.0-gam)*ebpts(:,kj+2,:,:);
            else
              ebpts(:,kj+1,:,:) = bet*ebpts(:,kj+1,:,:) + (1.0-bet)*ebpts(:,kj+2,:,:);
            end
          end
          i = i+1;
          j = j-1;
          kj = kj-1;
        end
        first = first-1;
        last = last+1;
      end
    end   % end of removing knot, u=V(a)
     
    if (a~=q+1)   % load the knot ua
      for i = 0:(qr-oldin-1)
         Vr(kind+1) = ua;
         kind = kind+1;
      end
    end
  
    for j = lbz:rbz   % load CPs into Qw
      Qw(:,cind+1,:,:) =  ebpts(:,j+1,:,:);
      cind = cind +1;
    end
  
    if (b<mv)   % setup for the next pass through loop
      for j = 0:in-1
        bpts(:,j+1,:,:) = Nextbpts(:,j+1,:,:);
      end
      for j = in:q
        bpts(:,j+1,:,:) = Pw(:,b-q+j,:,:);
      end
      a = b;
      b = b+1;
      ua = ub;
    else   % end knots
      for i = 0:qr
        Vr(kind+i+1) = ub;
      end
    end
  end
  % end of big loop through knot vector
elseif (tq==0)
  Vr=V;  qr=q;
end
%========================== DEGREE ELEVATE r ==============================
if (tr>0)
  clear ('Pw');
  Pw(:,:,:,:)=Qw(:,:,:,:);
  nu = length(Pw(:,1,1,1));
  nv = length(Pw(1,:,1,1));
  clear ('Qw');  
  clear ('bpts');
  clear ('Nextbpts');
  clear ('ebpts');
  clear ('bezalfs');
  
  mw = nw+r+1;
  rr = r+tr;
  % preallocate Wr
  kw = 1;
  for i=1:mw-1;  if(W(i)~=W(i+1));  kw = kw+1;   end;   end 
  Wr = zeros(1,mw+tr*kw);

  % compute Bezier degree elevation coefficients
  bezalfs(1,1) = 1.0;
  bezalfs(rr+1,r+1) = 1.0;

  for i = 1:fix(rr/2)   % bezalfs are symmetric to row rr/2
    inv = 1/nchoosek(rr,i);
    mpi = min(r,i);
    for j = max(0,i-tr):mpi   % alfas for degree elevation
      bezalfs(i+1,j+1) = inv*nchoosek(r,j)*nchoosek(tr,i-j);
      bezalfs(rr-i+1,r-j+1) = bezalfs(i+1,j+1);    % symmetry !
    end
  end

  mh = rr;
  kind = rr+1;
  in=-1;
  a=r+1;
  b=r+2;
  cind=1;
  ua = W(a);
  Qw(:,:,1,:) = Pw(:,:,1,:);

  % left end of Wr
  for i = 1:rr+1
    Wr(i) = ua;
  end

  % initialize first Bezier seg
  for i = 1:r+1
    bpts(:,:,i,:) = Pw(:,:,i,:);
  end

  % big loop through knot vector
  while (b<mw)
    i = b;
    while (b<mw && W(b)==W(b+1))
      b = b+1;      % make b the rightmost ocurrence of ub
    end
    mult = b-i+1;   % multiplicity of ub
    mh = mh+mult+tr;
    ub = W(b);
    oldin = in;       % in from last segment
    in = r-mult;     % ub to be inserted in times
  
    if (oldin>0);   lbz = fix((oldin+2)/2);
    else            lbz = 1;
    end
  
    if (in>0);      rbz = rr-fix((in+1)/2);
    else            rbz = rr;
    end
  
    % insert knots to get Bezier segment
    if (in>0)
      numer = ub - ua;
      alfs = zeros(1,r-mult);
      for k = r:-1:mult+1
       alfs(k-mult) = numer/(W(a+k)-ua);   % alfa for knot insertion
      end
      for j = 1:in   % in times knot insertion
        save = in-j+1;
        s = mult+j;
        for k = r+1:-1:s+1    % new CP due to knot insertion
          bpts(:,:,k,:) = alfs(k-s)*bpts(:,:,k,:)+(1-alfs(k-s))*bpts(:,:,k-1,:);
        end
        Nextbpts(:,:,save,:) = bpts(:,:,r+1,:);
      end
    end
    % end of inserting knots
  
    % degree elevate Bezier
    for i = lbz:rr
      % only points lbz,..,rr are used below
      ebpts(1:nu,1:nv,i+1,1:4) = 0;
      mpi = min(r,i);
      for j = max(0,i-tr):mpi   % new CP due to degree elevation
        ebpts(:,:,i+1,:) = ebpts(:,:,i+1,:) + bezalfs(i+1,j+1)*bpts(:,:,j+1,:);
      end
    end
    % end of degree elevating Bezier
  
    % knot removal ua oldin times
    if (oldin>1)
      first = kind-2;
      last = kind;
      den = ub-ua;
      bet = (ub-Wr(kind))/den;
      for tr = 1:(oldin-1)
        i = first;
        j = last;
        kj = j-kind+1;
        while ((j-i)>tr)
          % loop and compute the new CP for one removal step
          if (i<cind)
            alf = (ub-Wr(i+1))/(ua-Wr(i+1));
            Qw(:,:,i+1,:) = (alf*Qw(:,:,i+1,:)+(1.0-alf)*Qw(:,:,i,:));
          end
          if (j>=lbz)
            if ((j-tr)<=(kind-rr+oldin))
              gam = (ub-Wr(j-tr+1))/den;
              ebpts(:,:,kj+1,:) = gam*ebpts(:,:,kj+1,:) + (1.0-gam)*ebpts(:,:,kj+2,:);
            else
              ebpts(:,:,kj+1,:) = bet*ebpts(:,:,kj+1,:) + (1.0-bet)*ebpts(:,:,kj+2,:);
            end
          end
          i = i+1;
          j = j-1;
          kj = kj-1;
        end
        first = first-1;
        last = last+1;
      end
    end   % end of removing knot, u=W(a)
     
    if (a~=r+1)   % load the knot ua
      for i = 0:(rr-oldin-1)
         Wr(kind+1) = ua;
         kind = kind+1;
      end
    end
  
    for j = lbz:rbz   % load CPs into Qw
      Qw(:,:,cind+1,:) =  ebpts(:,:,j+1,:);
      cind = cind +1;
    end
  
    if (b<mw)   % setup for the next pass through loop
      for j = 0:in-1
        bpts(:,:,j+1,:) = Nextbpts(:,:,j+1,:);
      end
      for j = in:r
        bpts(:,:,j+1,:) = Pw(:,:,b-r+j,:);
      end
      a = b;
      b = b+1;
      ua = ub;
    else   % end knots
      for i = 0:rr
        Wr(kind+i+1) = ub;
      end
    end
  end
  % end of big loop through knot vector
elseif (tr==0)
  Wr=W;  rr=r;
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