% 
% J. Kiendl
clear;

deg = 1;
ref = 64;
slen = 100;
L=100; t=L/slen; w=1;
% material
E = 1e3;
nue = 0.0;

% Gauss Points
ngauss=[2 2 2];

p=1; q=1; r=1;
U=[0 0 1 1]; V=[0 0 1 1]; W=[0 0 1 1];
clear CP;
CP(:,:,1,1) = [0 0; 1 1]*L;
CP(:,:,1,2) = [0 1; 0 1]*w;
CP(:,:,1,3) = [0 0; 0 0];
CP(:,:,1,4) = [1 1; 1 1];
CP(:,:,2,:) = CP(:,:,1,:);
CP(:,:,2,3) = CP(:,:,1,3)+t;
% plotNURBS_solid(p,q,r,U,V,W,CP)

[p,q,r,U,V,W,CP] = degree_elevate_solid(p,q,r,U,V,W,CP,deg-p,deg-q,deg-r);
R = refinement_vec(U,ref);
[CP,U,V,W] = knot_refine_solid(p,q,r,U,V,W,CP,R,[],[]);
% plotNURBS_solid(p,q,r,U,V,W,CP)
% fixed supports
rb =[];
us=[0 0]; vs=[0 1]; ws=[0 1];  dirs=1;
rb = supports(rb,us,vs,ws,dirs,CP);
us=[0 0]; vs=[0 1]; ws=[0 1];  dirs=2;
rb = supports(rb,us,vs,ws,dirs,CP);
us=[0 0]; vs=[0 1]; ws=[0 1];  dirs=3;
rb = supports(rb,us,vs,ws,dirs,CP);

% load vector
F = 0.1*t^(1/2);
%F = E*w/4/L^3*t^3;  
fl=[];
ub=1; vb=[0 1]; wb=[0 1];  dirf=3;
fl = load_area(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,F,dirf,0);
% fl2 = load_area(fl2,ub,vb,wb,p2,q2,r2,U2,V2,W2,CP2,ngauss,-fq,dirf,proj);

% compute displacement vector d and load vector including supports fs
  [dc,fs] = solve(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);
% [da,fs] = solve_ANS(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);
  CPdc = CPresult(CP,dc);
%   CPda = CPresult(CP,da);

en = 0.5*fs'*dc;
disp(en);

% visualize
  plot2in1(p,q,r,U,V,W,CP,CPdc,rb,fl);
% plot_strain_stress(p,q,r,U,V,W,CP,CPd,rb,fl,E,nue,d,2,2)

