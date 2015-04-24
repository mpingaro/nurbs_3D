% 
clear;
p=2;
q=2;
r=2;
U=[0 0 0 1 1 1];
V=[0 0 0 1 1 1];
W=[0 0 0 1 1 1];
CP(:,:,1,1) = [0 0 0; 1 1 1; 2 2 2];
CP(:,:,2,1) = [3 1 -1; 3 1 -1; 3 1 -1];
CP(:,:,3,1) = [2 2 2; 1 1 1; 0 0 0];
CP(:,:,1,2) = [0 1 2; 0 1 2; 0 1 2];
CP(:,:,2,2) = [-1 -1 -1; 1 1 1; 3 3 3];
CP(:,:,3,2) = [2 1 0; 2 1 0; 2 1 0];
CP(:,:,1,3) = [0 0 0; 0 0 0; 0 0 0];
CP(:,:,2,3) = [1 1 1; 1 1 1; 1 1 1]*5;
CP(:,:,3,3) = [1 1 1; 1 1 1; 1 1 1]*10;
CP(:,:,1,4) = [1 1 1; 1 1 1; 1 1 1];
CP(:,:,2,4) = [1 1 1; 1 1 1; 1 1 1];
CP(:,:,3,4) = [1 1 1; 1 1 1; 1 1 1];
plotNURBS_solid(p,q,r,U,V,W,CP)
% material
E = 1e3;
nue = 0.0;

% Gauss Points
ngauss=[3 3 3];

Ru=[0.5]; Rv=[0.5]; Rw=[0.2 0.4 0.6 0.8];
[CP,U,V,W] = knot_refine_solid(p,q,r,U,V,W,CP,Ru,Rv,Rw);
plotNURBS_solid(p,q,r,U,V,W,CP)
% fixed supports
rb =[];
us=[0 0]; vs=[0 1]; ws=[0 0];  dirs=1;
rb = supports(rb,us,vs,ws,dirs,CP);
us=[0 1]; vs=[0 1]; ws=[0 0];  dirs=2;
rb = supports(rb,us,vs,ws,dirs,CP);
us=[0 1]; vs=[0 1]; ws=[0 0];  dirs=3;
rb = supports(rb,us,vs,ws,dirs,CP);

% load vector
F = -1.5;  fl=[];
ub=[0 1]; vb=[0 1]; wb=1;  dirf=1;
fl = load_area(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,F,dirf,0);
% fl2 = load_area(fl2,ub,vb,wb,p2,q2,r2,U2,V2,W2,CP2,ngauss,-fq,dirf,proj);

% compute displacement vector d and load vector including supports fs
[d,fs] = solve(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);
CPd = CPresult(CP,d);

% visualize
plot2in1(p,q,r,U,V,W,CP,CPd,rb,fl)
plot_strain_stress(p,q,r,U,V,W,CP,CPd,rb,fl,E,nue,d,2,2)