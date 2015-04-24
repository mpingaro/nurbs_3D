% Example from ANS solid-shell paper
% 2. cantilever beam with 8 Element and varying thickness
clear;

SL = 5
deg = 3;

p=2;
q=2;
r=2;
U=[0 0 0 1 1 1];
V=[0 0 0 0.125 0.25 0.375 0.5 0.625 0.75 0.875 1 1 1];
W=[0 0 0 1 1 1];

fold = '/Users/kiendljosef/Desktop/Forschung/solid_shell/Input_Files_edited/';
name = strcat('BernoulliStrBeam_SL',int2str(SL),'(2)_p2q2w2');
CPm = textread(strcat(fold,name));
CP = CP_matrix(p,q,r,U,V,W,CPm);

[p,q,r,U,V,W,CP] = degree_elevate_solid(p,q,r,U,V,W,CP,deg-p,deg-q,deg-r);

% material
E = 1e3;
nue = 0.0;
% Gauss Points
ngauss=[p+1 q+1 r+1];

% fixed supports
rb =[];
us=[0 1]; vs=[0 0]; ws=[0 1];  dirs=1;
rb = supports(rb,us,vs,ws,dirs,CP);
us=[0 1]; vs=[0 0]; ws=[0 1];  dirs=2;
rb = supports(rb,us,vs,ws,dirs,CP);
us=[0 1]; vs=[0 0]; ws=[0 1];  dirs=3;
rb = supports(rb,us,vs,ws,dirs,CP);

% load vector
f = 0.1*(CPm(3,3)-CP(1,3))^(1/2);  fl=[];
% ub=[0 1]; vb=[0 1]; wb=[0 1];  dirf=3;
% fl = load_body(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,f,dirf);
ub=[0 1]; vb=1; wb=[0 1];  dirf=3;  proj=0;
fl = load_area(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,f,dirf,proj);
% ub=[0 1]; vb=0; wb=1;  dirf=3;  proj=0;
% fl = load_line(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,f,dirf,proj);

% compute displacement vector d and load vector including supports fs
[d,fs] = solve(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);
CPd = CPresult(CP,d);

en = 0.5*fl'*d

% name2 = strcat('_deg',int2str(deg),'ref',int2str(ref));
% name3 = strcat(fold,name,name2,'.mat');
% save(name3,'en');

% visualize
% plot2in1(p,q,r,U,V,W,CP,CPd,rb,fl)