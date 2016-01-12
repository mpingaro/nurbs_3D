function plot_strain_stress(p,q,r,U,V,W,CP1,CP2,rb,f,E,nue,d,ss,comp)
% plots two windows. The first one shows the undeformed geometry with
% loads and supports. The second one shows the strains or stresses on the
% deformed geometry.
% Parameters:
% CP1: undeformed control points
% CP2: deformed control points
% ss: 1=plot strain;   2=plot stress
% comp: component of strain/stress:  1-6: xx,yy,zz,xy,yz,zx
%                                    7-9: principal strains/stresses
%                                     10: von Mises stress

mu = length(U);
mv = length(V);
mw = length(W);
nu = length(CP2(:,1,1,1));
nv = length(CP2(1,:,1,1));
nw = length(CP2(1,1,:,1));

D = E/((1+nue)*(1-2*nue)).*[1-nue nue nue 0 0 0; nue 1-nue nue 0 0 0; nue nue 1-nue 0 0 0;...
    0 0 0 (1-2*nue)/2 0 0; 0 0 0 0 (1-2*nue)/2 0; 0 0 0 0 0 (1-2*nue)/2];
%=================================================================
% FIRST WINDOW: UNDEFORMED STRUCTURE
% 1. original geometry, supports and load arrows
[xs1,ys1,zs1] = create_supports(CP1,rb);
[xf1,yf1,zf1] = create_arrows(CP1,f);

% FIRST WINDOW: UNDEFORMED
subplot(2,1,1);

% geometry
plotNURBS_solid(p,q,r,U,V,W,CP1);
hold;

% supports
for k =1:length(xs1(:,1))
  plot3(xs1(k,:),ys1(k,:),zs1(k,:),'Linewidth',2,'Color','black');
end

% load arrows
for k =1:length(xf1(:,1))
  plot3(xf1(k,:),yf1(k,:),zf1(k,:),'Linewidth',5);
  plot3(xf1(k,1),yf1(k,1),zf1(k,1),'Marker','d','MarkerFaceColor','blue','MarkerSize',10);
end

axis equal;
camlight left;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('z','FontSize',14);
title ('original geometry');
hold off;

%=================================================================
% SECOND WINDOW: STRAIN/STRESS ON DEFORMED STRUCTURE
subplot(2,1,2);
% element displacement vectors
d_el = zeros(mu-p-1,mv-q-1,mw-r-1,3*(p+1)*(q+1)*(r+1));
for k = (r+1):(mw-r-1)
  for j = (q+1):(mv-q-1)
    for i = (p+1):(mu-p-1)
      l=1; 
      for c = k-r-1:k-1
        for b = j-q-1:j-1 
          for a = i-p-1:i-1
              d_el(i,j,k,l)   = d(3*(c*nu*nv+b*nu+a)+1);
              d_el(i,j,k,l+1) = d(3*(c*nu*nv+b*nu+a)+2);
              d_el(i,j,k,l+2) = d(3*(c*nu*nv+b*nu+a)+3);
              l=l+3;
          end
        end
      end
    end
  end
end
% W =======================================================================
grid = 20;
tol=10e-10;
w=W(1);
for wi = 1:nw-1:nw
  k = findspan(w,W,nw);
  sv=(V(mv)-V(1))/grid;
  v=V(1);
  b=1;
  while v <= V(mv)+tol
    j = findspan(v,V,nv);
    su=(U(mu)-U(1))/grid;
    u=U(1);
    a=1;
    while u <= U(mu)+tol
      i = findspan(u,U,nu);
      P(a,b,1:3) = get_point_solid(u,v,w,i,j,k,p,q,r,U,V,W,CP2);
      d_actual(:,1) = d_el(i,j,k,:);
      eps_actual = strain(p,q,r,i,j,k,u,v,w,U,V,W,CP1,d_actual);
      % principal strains from eigenvalues of strain tensor
      EPS = [    eps_actual(1) 0.5*eps_actual(4) 0.5*eps_actual(6);
             0.5*eps_actual(4)     eps_actual(2) 0.5*eps_actual(5);
             0.5*eps_actual(6) 0.5*eps_actual(5)     eps_actual(3)];
      e = eig(EPS);   e = sort(e,'descend');
      eps(a,b,1:6) = eps_actual;
      eps(a,b,7:9) = e(1:3);
      % stresses
      sig(a,b,1:6) = D*eps_actual;
      % principal stresses from principal strains
      % e = principal strains vector with shear strains equal 0.
      e(4:6) = 0;
      s = D*e;
      sig(a,b,7:9) = s(1:3);
      % von Mises stress from principal stresses
      sig(a,b,10) = sqrt(s(1)^2+s(2)^2+s(3)^2-s(1)*s(2)-s(2)*s(3)-s(3)*s(1));
      a=a+1;
      u=u+su;
    end
    b=b+1;
    v=v+sv;
  end
  if     (ss==1); surf(P(:,:,1),P(:,:,2),P(:,:,3),eps(:,:,comp));
  elseif (ss==2); surf(P(:,:,1),P(:,:,2),P(:,:,3),sig(:,:,comp));   end
  shading interp;
  hold on;
  CP_surf(:,:,:) = CP2(:,:,wi,:);
  create_el_edges(p,q,U,V,CP_surf);
  w=W(mw);
end
clear P; clear eps; clear sig; clear CP_surf;
% V =======================================================================
v=V(1);
for vi = 1:nv-1:nv
  j = findspan(v,V,nv);
  sw=(W(mw)-W(1))/grid;
  w=W(1);
  b=1;
  while w <= W(mw)+tol
    k = findspan(w,W,nw);
    su=(U(mu)-U(1))/grid;
    u=U(1);
    a=1;
    while u <= U(mu)+tol
      i = findspan(u,U,nu);
      P(a,b,1:3) = get_point_solid(u,v,w,i,j,k,p,q,r,U,V,W,CP2);
      d_actual(:,1) = d_el(i,j,k,:);
      eps_actual = strain(p,q,r,i,j,k,u,v,w,U,V,W,CP1,d_actual);
      % principal strains from eigenvalues of strain tensor
      EPS = [    eps_actual(1) 0.5*eps_actual(4) 0.5*eps_actual(6);
             0.5*eps_actual(4)     eps_actual(2) 0.5*eps_actual(5);
             0.5*eps_actual(6) 0.5*eps_actual(5)     eps_actual(3)];
      e = eig(EPS);   e = sort(e,'descend');
      eps(a,b,1:6) = eps_actual;
      eps(a,b,7:9) = e(1:3);
      % stresses
      sig(a,b,1:6) = D*eps_actual;
      % principal stresses from principal strains
      % e = principal strains vector with shear strains equal 0.
      e(4:6) = 0;
      s = D*e;
      sig(a,b,7:9) = s(1:3);
      % von Mises stress
      sig(a,b,10) = sqrt(s(1)^2+s(2)^2+s(3)^2-s(1)*s(2)-s(2)*s(3)-s(3)*s(1));
      a=a+1;
      u=u+su;
    end
    b=b+1;
    w=w+sw;
  end
  if     (ss==1); surf(P(:,:,1),P(:,:,2),P(:,:,3),eps(:,:,comp));
  elseif (ss==2); surf(P(:,:,1),P(:,:,2),P(:,:,3),sig(:,:,comp));   end
  shading interp;
  hold on;
  CP_surf(:,:,:) = CP2(:,vi,:,:);
  create_el_edges(p,r,U,W,CP_surf);
  v=V(mv);
end
clear P; clear eps; clear sig; clear CP_surf;
% U =======================================================================
u=U(1);
for ui = 1:nu-1:nu
  i = findspan(u,U,nu);
  sw=(W(mw)-W(1))/grid;
  w=W(1);
  b=1;
  while w <= W(mw)+tol
    k = findspan(w,W,nw);
    sv=(V(mv)-V(1))/grid;
    v=V(1);
    a=1;
    while v <= V(mv)+tol
      j = findspan(v,V,nv);
      P(a,b,1:3) = get_point_solid(u,v,w,i,j,k,p,q,r,U,V,W,CP2);
      d_actual(:,1) = d_el(i,j,k,:);
      eps_actual(:,1) = strain(p,q,r,i,j,k,u,v,w,U,V,W,CP1,d_actual);
      % principal strains from eigenvalues of strain tensor
      EPS = [    eps_actual(1) 0.5*eps_actual(4) 0.5*eps_actual(6);
             0.5*eps_actual(4)     eps_actual(2) 0.5*eps_actual(5);
             0.5*eps_actual(6) 0.5*eps_actual(5)     eps_actual(3)];
      e = eig(EPS);   e = sort(e,'descend');
      eps(a,b,1:6) = eps_actual;
      eps(a,b,7:9) = e(1:3);
      % stresses
      sig(a,b,1:6) = D*eps_actual;
      % principal stresses from principal strains
      % e = principal strains vector with shear strains equal 0.
      e(4:6) = 0;
      s = D*e;
      sig(a,b,7:9) = s(1:3);
      % von Mises stress
      sig(a,b,10) = sqrt(s(1)^2+s(2)^2+s(3)^2-s(1)*s(2)-s(2)*s(3)-s(3)*s(1));
      a=a+1;
      v=v+sv;
    end
    b=b+1;
    w=w+sw;
  end
  if     (ss==1); surf(P(:,:,1),P(:,:,2),P(:,:,3),eps(:,:,comp));
  elseif (ss==2); surf(P(:,:,1),P(:,:,2),P(:,:,3),sig(:,:,comp));   end
  shading interp;
  hold on;
  CP_surf(:,:,:) = CP2(ui,:,:,:);
  create_el_edges(q,r,V,W,CP_surf);
  u=U(mu);
end
clear P; clear eps; clear sig; clear CP_surf;
% =========================================================================
colormap('default');
% invert default colormap => red=negativ, blue=positive
COL=colormap;
invCOL(:,1)=COL(:,3);
invCOL(:,2)=COL(:,2);
invCOL(:,3)=COL(:,1);
colormap(invCOL);
% make colormap symmetric
colim = caxis;
caxis([-max(abs(colim)) max(abs(colim))]);
% ------
colorbar;
axis equal;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('z','FontSize',14);
if     (ss==1 && comp==1);  title ('strain epsilon xx');
elseif (ss==1 && comp==2);  title ('strain epsilon yy');
elseif (ss==1 && comp==3);  title ('strain epsilon zz');
elseif (ss==1 && comp==4);  title ('strain epsilon xy');
elseif (ss==1 && comp==5);  title ('strain epsilon yz');
elseif (ss==1 && comp==6);  title ('strain epsilon zx');
elseif (ss==1 && comp==7);  title ('principle strain e1');
elseif (ss==1 && comp==8);  title ('principle strain e2');
elseif (ss==1 && comp==9);  title ('principle strain e3');
elseif (ss==2 && comp==1);  title ('stress sigma xx');
elseif (ss==2 && comp==2);  title ('stress sigma yy');
elseif (ss==2 && comp==3);  title ('stress sigma zz');
elseif (ss==2 && comp==4);  title ('stress sigma xy');
elseif (ss==2 && comp==5);  title ('stress sigma yz');
elseif (ss==2 && comp==6);  title ('stress sigma zx');
elseif (ss==2 && comp==7);  title ('principle stress s1');
elseif (ss==2 && comp==8);  title ('principle stress s2');
elseif (ss==2 && comp==9);  title ('principle stress s3');
elseif (ss==2 && comp==10); title ('von Mises stress');
end
hold off;
