% NURBS parameters
clear;
p=2;
q=1;
r=1;
U=[0 0 0 1 1 1];
V=[0 0 1 1];
W=[0 0 1 1];
load 'CPstiffener';
CP1=CP3d(:,1:2,:,:);
CP2=CP3d(:,2:3,:,:);

% plotNURBS_solid(p,q,r,U,V,W,CP1); hold on;
% plotNURBS_solid(p,q,r,U,V,W,CP2);

% material
E = 6.825e7;
nue = 0.3;

deg=3;
refu=4;
refv=1;
refw=2;

% Gauss Points
ngauss=[deg+2 deg+2 deg+2];


[p1,q1,r1,U1,V1,W1,CP1] = degree_elevate_solid(p,q,r,U,V,W,CP1,deg-p,deg-q,deg-r);
[p2,q2,r2,U2,V2,W2,CP2] = degree_elevate_solid(p,q,r,U,V,W,CP2,deg-p,deg-q,deg-r);
Ru=[]; Rv=[]; Rw=[];
for i = 1:refu-1;  Ru(i) = i/refu;  end
for i = 1:refv-1;  Rv(i) = i/refv;  end
for i = 1:refw-1;  Rw(i) = i/refw;  end
[CP1,U1,V1,W1] = knot_refine_solid(p1,q1,r1,U1,V1,W1,CP1,Ru,Rv,Rw);
[CP2,U2,V2,W2] = knot_refine_solid(p2,q2,r2,U2,V2,W2,CP2,Ru,Rv,Rw);

nu1 = length(CP1(:,1,1,1)); nv1 = length(CP1(1,:,1,1)); nw1 = length(CP1(1,1,:,1));
ndof1=3*nu1*nv1*nw1;
nu2 = length(CP2(:,1,1,1)); nv2 = length(CP2(1,:,1,1)); nw2 = length(CP2(1,1,:,1));
ndof2=3*nu2*nv2*nw2;
ndof12=ndof1+ndof2;

% plotNURBS_solid(p1,q1,r1,U1,V1,W1,CP1); hold on;
% plotNURBS_solid(p2,q2,r2,U2,V2,W2,CP2);

K1 = stiff_mat(p1,U1,q1,V1,r1,W1,CP1,E,nue,ngauss);
K2 = stiff_mat(p2,U2,q2,V2,r2,W2,CP2,E,nue,ngauss);

% fixed supports
rb1 = supports([],[0 1],[0 1],[1 1],3,CP1);
rb1 = supports(rb1,[0 0],[0 1],[0 1],2,CP1);
rb1 = supports(rb1,[1 1],[0 1],[0 1],1,CP1);
rb2 = rb1;

% load
fg=-5; fq=100;
fl1=[];
ub=[0 1]; vb=[0 1]; wb=[0 1];  dirf=3;
fl1 = load_body(fl1,ub,vb,wb,p1,q1,r1,U1,V1,W1,CP1,ngauss,fg,dirf);
fl2=[];
% ub=[0 1]; vb=[0 1]; wb=0;  dirf=7;  proj=0;
% fl2 = load_area(fl2,ub,vb,wb,p2,q2,r2,U2,V2,W2,CP2,ngauss,fq,dirf,proj);
ub=[0 1]; vb=1; wb=[0 1];  dirf=7;  proj=0;
fl2 = load_area(fl2,ub,vb,wb,p2,q2,r2,U2,V2,W2,CP2,ngauss,-fq,dirf,proj);
ub=[0 1]; vb=[0 1]; wb=[0 1];  dirf=3;
fl2 = load_body(fl2,ub,vb,wb,p2,q2,r2,U2,V2,W2,CP2,ngauss,fg,dirf);

save Kstiffener_p3_r4 CP1 CP2 K1 K2 U1 U2 V1 V2 W1 W2 p1 p2 q1 q2 r1 r2 ...
     rb1 rb2 fl1 fl2 ndof1 ndof2 nu1 nv1 nw1 nu2 nv2 nw2

load 'Kstiffener_d_p2'

% scal=500;
d1=scal*d(1:ndof1);         CP1d = CPresult(CP1,d1);
d2=scal*d(ndof1+1:ndof12);  CP2d = CPresult(CP2,d2);

plotNURBS_solid(p1,q1,r1,U1,V1,W1,CP1d); hold on;
plotNURBS_solid(p2,q2,r2,U2,V2,W2,CP2d);