function plotNURBS_solid(p,q,r,U,V,W,CP)
% plots the six faces of a solid
mu=length(U); mv=length(V); mw=length(W);
nu=length(CP(:,1,1,1));
nv=length(CP(1,:,1,1));
nw=length(CP(1,1,:,1));

check_input(p,mu,nu,q,mv,nv,r,mw,nw)
plot3(0,0,0);  hold on;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('z','FontSize',14);
% W ---------------------
CP1=zeros(nu,nv,4);
CP1(:,:,:)=CP(:,:,1,:);
plotNURBS_surf(p,q,U,V,CP1);
create_el_edges(p,q,U,V,CP1);
CP1(:,:,:)=CP(:,:,nw,:);
plotNURBS_surf(p,q,U,V,CP1);
create_el_edges(p,q,U,V,CP1);
% V ---------------------
CP2=zeros(nu,nw,4);
CP2(:,:,:)=CP(:,1,:,:);
plotNURBS_surf(p,r,U,W,CP2);
create_el_edges(p,r,U,W,CP2);
CP2(:,:,:)=CP(:,nv,:,:);
plotNURBS_surf(p,r,U,W,CP2);
create_el_edges(p,r,U,W,CP2);
% U ---------------------
CP3=zeros(nv,nw,4);
CP3(:,:,:)=CP(1,:,:,:);
plotNURBS_surf(q,r,V,W,CP3);
create_el_edges(q,r,V,W,CP3);
CP3(:,:,:)=CP(nu,:,:,:);
plotNURBS_surf(q,r,V,W,CP3);
create_el_edges(q,r,V,W,CP3);
% -----------------------
% create_conpoints(CP)
% create_conpolygon(CP)
% -----------------------
axis equal;
camlight left;
% lighting phong;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('z','FontSize',14);
hold off;