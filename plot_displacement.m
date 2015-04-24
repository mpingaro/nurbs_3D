function plot_displacement(p,q,r,U,V,W,CP1,CP2,rb,f,dir)
% plots displacement contours

% Parameters:
% CP1: undeformed control points
% CP2: deformed control points
% dir: direction:  1-x  2-y  3-z  4-radial (for cylinder example, not general)

mu = length(U);
mv = length(V);
mw = length(W);
nu = length(CP2(:,1,1,1));
nv = length(CP2(1,:,1,1));
nw = length(CP2(1,1,:,1));

% structure
% W =======================================================================
grid = 100;
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
      P1(a,b,1:3) = get_point_solid(u,v,w,i,j,k,p,q,r,U,V,W,CP1);
      P2(a,b,1:3) = get_point_solid(u,v,w,i,j,k,p,q,r,U,V,W,CP2);
      d(a,b,1:3)  = P2(a,b,1:3) - P1(a,b,1:3);
      a=a+1;
      u=u+su;
    end
    b=b+1;
    v=v+sv;
  end
  % for cylinder example: refer x to the new center
  dxm(:,:) = d(:,:,1) - (CP2(1,1,1,1)+CP2((nu+1)/2,1,1,1))/2;
  d(:,:,4) = sqrt(dxm(:,:).^2 + d(:,:,2).^2);
  surf(P2(:,:,1),P2(:,:,2),P2(:,:,3),d(:,:,dir));
  shading interp;
  hold on;
  CP_surf(:,:,:) = CP2(:,:,wi,:);
  create_el_edges(p,q,U,V,CP_surf);
  w=W(mw);
end
clear P1; clear P2; clear CP_surf;
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
      P1(a,b,1:3) = get_point_solid(u,v,w,i,j,k,p,q,r,U,V,W,CP1);
      P2(a,b,1:3) = get_point_solid(u,v,w,i,j,k,p,q,r,U,V,W,CP2);
      d(a,b,1:3)  = P2(a,b,1:3) - P1(a,b,1:3);
      a=a+1;
      u=u+su;
    end
    b=b+1;
    w=w+sw;
  end
  % for cylinder example: refer x to the new center
  dxm(:,:) = d(:,:,1) - (CP2(1,1,1,1)+CP2((nu+1)/2,1,1,1))/2;
  d(:,:,4) = sqrt(dxm(:,:).^2 + d(:,:,2).^2);
  surf(P2(:,:,1),P2(:,:,2),P2(:,:,3),d(:,:,dir));
  shading interp;
  hold on;
  CP_surf(:,:,:) = CP2(:,vi,:,:);
  create_el_edges(p,r,U,W,CP_surf);
  v=V(mv);
end
clear P1; clear P2; clear CP_surf;
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
      P1(a,b,1:3) = get_point_solid(u,v,w,i,j,k,p,q,r,U,V,W,CP1);
      P2(a,b,1:3) = get_point_solid(u,v,w,i,j,k,p,q,r,U,V,W,CP2);
      d(a,b,1:3)  = P2(a,b,1:3) - P1(a,b,1:3);
      a=a+1;
      v=v+sv;
    end
    b=b+1;
    w=w+sw;
  end
  % for cylinder example: refer x to the new center
  dxm(:,:) = d(:,:,1) - (CP2(1,1,1,1)+CP2((nu+1)/2,1,1,1))/2;
  d(:,:,4) = sqrt(dxm(:,:).^2 + d(:,:,2).^2);
  surf(P2(:,:,1),P2(:,:,2),P2(:,:,3),d(:,:,dir));
  shading interp;
  hold on;
  CP_surf(:,:,:) = CP2(ui,:,:,:);
  create_el_edges(q,r,V,W,CP_surf);
  u=U(mu);
end
clear P1; clear P2; clear CP_surf;
% =========================================================================
colormap('default');
% % invert default colormap => red=negativ, blue=positive
% COL=colormap;
% invCOL(:,1)=COL(:,3);
% invCOL(:,2)=COL(:,2);
% invCOL(:,3)=COL(:,1);
% colormap(invCOL);
% % make colormap symmetric
% colim = caxis;
% caxis([-max(abs(colim)) max(abs(colim))]);
% ------
colorbar;
axis equal;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('z','FontSize',14);
if     (dir==1);  title ('displacement in x','FontSize',20);
elseif (dir==2);  title ('displacement in y','FontSize',20);
elseif (dir==3);  title ('displacement in z','FontSize',20);
elseif (dir==4);  title ('radial displacement','FontSize',20);
end
%=================================================================
[xs,ys,zs] = create_supports(CP2,rb);
[xf,yf,zf] = create_arrows(CP2,f);

% supports
% for k =1:length(xs(:,1))
%   plot3(xs(k,:),ys(k,:),zs(k,:),'Linewidth',2,'Color','black');
%   hold on;
% end

% load arrows
% for k =1:length(xf(:,1))
%   plot3(xf(k,:),yf(k,:),zf(k,:),'Linewidth',5);
%   plot3(xf(k,1),yf(k,1),zf(k,1),'Marker','d','MarkerFaceColor','blue','MarkerSize',10);
% end
%=================================================================
hold off;