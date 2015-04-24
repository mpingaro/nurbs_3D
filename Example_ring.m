% Cylinder
clear;
p=2;
q=1;
r=1;

U = [0 0 0 1 1 2 2 3 3 4 4 4];          
V = [0 0 1 1];  
W = [0 0 1 1];

CP(:,:,1,1) = [10 12; 10 12;  0  0; -10 -12; -10 -12; -10 -12;   0   0;  10  12;  10 12];
CP(:,:,2,1) = [10 12; 10 12;  0  0; -10 -12; -10 -12; -10 -12;   0   0;  10  12;  10 12];

CP(:,:,1,2) = [ 0  0; 10 12; 10 12;  10  12;   0   0; -10 -12; -10 -12; -10 -12;  0   0];
CP(:,:,2,2) = [ 0  0; 10 12; 10 12;  10  12;   0   0; -10 -12; -10 -12; -10 -12;  0   0];

CP(:,:,1,3) =  [0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0];
CP(:,:,2,3) = -[1 1; 1 1; 1 1; 1 1; 1 1; 1 1; 1 1; 1 1; 1 1];

w2=sqrt(0.5);
CP(:,:,1,4) = [1 1; w2 w2; 1 1; w2 w2; 1 1; w2 w2; 1 1; w2 w2; 1 1];
CP(:,:,2,4) = [1 1; w2 w2; 1 1; w2 w2; 1 1; w2 w2; 1 1; w2 w2; 1 1];

figure(1);
plotNURBS_solid(p,q,r,U,V,W,CP);
axis off;

[p,q,r,U,V,W,CP] = degree_elevate_solid(p,q,r,U,V,W,CP,1,2,0);
% Ru=[0.5 1.5 2.5 3.5]; Rv=[]; Rw=[];
Ru=[]; Rv=[]; Rw=[];
[CP,U,V,W] = knot_refine_solid(p,q,r,U,V,W,CP,Ru,Rv,Rw);

% plotNURBS_solid(p,q,r,U,V,W,CP)
% material
E = 1e3;  % E = 7.5e7;
nue = 0.0;
% Gauss Points
ngauss=[p+1 q+1 r+1];

% closed faces
closed = [1 0 0];
ms = master_slaves(closed,CP);

% fixed supports
rb =[];
us=[0 0]; vs=[0 0]; ws=[0 0];  dirs=1;
rb = supports(rb,us,vs,ws,dirs,CP);
us=[0 0]; vs=[0 1]; ws=[0 0];  dirs=2;
rb = supports(rb,us,vs,ws,dirs,CP);
us=[0 1]; vs=[0 1]; ws=[0 0];  dirs=3;
rb = supports(rb,us,vs,ws,dirs,CP);

% load vector
f = 1;  fl=[];
ub=[0 4]; vb=0; wb=[0 1];  dirf=7;   proj=0;
fl = load_area(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,f,dirf,proj);

% compute displacement vector d and load vector including supports fs
% [d,fs] = solve(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);
[d,fs] = solve_closed(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb,ms);
CPd = CPresult(CP,d);

% visualize
% plot2in1(p,q,r,U,V,W,CP,CPd,rb,fl);
% plot_strain_stress(p,q,r,U,V,W,CP,CPd,E,nue,d,2,9);
figure(2);
plot_strain_stress_1in1(p,q,r,U,V,W,CP,CPd,rb,fl,E,nue,d,2,9);