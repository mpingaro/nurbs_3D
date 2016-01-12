% Cylinder
clear;
clc;
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

% figure(1);
% plotNURBS_solid(p,q,r,U,V,W,CP);
% axis off;

[p,q,r,U,V,W,CP] = degree_elevate_solid(p,q,r,U,V,W,CP,1,1,1);
% Ru=[0.5 1.5 2.5 3.5]; Rv=[]; Rw=[];
% Ru=[]; Rv=[]; Rw=[];
Ru = refinement_vec(U,2);
Rv = refinement_vec(V,1);
Rw = refinement_vec(W,1);
[CP,U,V,W] = knot_refine_solid(p,q,r,U,V,W,CP,Ru,Rv,Rw);

mu = length(U);
mv = length(V);
mw = length(W);
nu = length(CP(:,1,1,1));
nv = length(CP(1,:,1,1));
nw = length(CP(1,1,:,1));
nn = nu*nv*nw;
ndof = 3*nn;

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
cp_rb1 = [1:(nu-1)/4+1];
cp_rb2 = cp_rb1;
for j=1:nv-1; cp_rb2 = [cp_rb2 j*nu+cp_rb1]; end
cp_rb = cp_rb2;
for k=1:nw-1; cp_rb = [cp_rb k*nu*nv+cp_rb2]; end
rb = sort([3*cp_rb-2 3*cp_rb-1 3*cp_rb]);
% rb =[];
% us=[0 1]; vs=[0 1]; ws=[0 0];  dirs=3;
% rb = supports(rb,us,vs,ws,dirs,CP);
% rb = rb(1:round(length(rb)/2));
% us=[0 0]; vs=[0 0]; ws=[0 0];  dirs=1;
% rb = supports(rb,us,vs,ws,dirs,CP);
% us=[0 0]; vs=[0 1]; ws=[0 0];  dirs=2;
% rb = supports(rb,us,vs,ws,dirs,CP);

% load vector
fl=[];
% f = -0.03;  
% ub=[2 3]; vb=[0 1]; wb=0;  dirf=3;   proj=0;
% fl = load_area(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,f,dirf,proj);
f = -0.5;
ub=2; vb=[0 1]; wb=0;  dirf=3;   proj=0;
fl = load_line(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,f,dirf,proj);
f = 0.5;
ub=3; vb=[0 1]; wb=0;  dirf=3;   proj=0;
fl = load_line(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,f,dirf,proj);
% fl = zeros(ndof,1);
% fl()

% plot1in1(p,q,r,U,V,W,CP,rb,fl);

% compute displacement vector d and load vector including supports fs
% [d,fs] = solve(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);
[d,fs] = solve_closed(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb,ms);
CPd = CPresult(CP,d);

% visualize
% plot1in1(p,q,r,U,V,W,CPd,rb,fl);
% plot2in1(p,q,r,U,V,W,CP,CPd,rb,fl);
 plot_displacement(p,q,r,U,V,W,CP,CPd,rb,fl,3); axis off;
%plot_strain_stress(p,q,r,U,V,W,CP,CPd,E,nue,d,2,9);
% figure(2);
% plot_strain_stress_1in1(p,q,r,U,V,W,CP,CPd,rb,fl,E,nue,d,2,7);