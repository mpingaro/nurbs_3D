% Example from ANS solid-shell paper
% 4. Scordelis-Lo
clear;

% deg = 2;
refs = [2 4 8 16];

for deg = 2:2
for iref = 1:length(refs)
  ref = refs(iref);

p=2;
q=2;
r=2;
U=[0 0 0 0.5 1 1 1];
V=[0 0 0 0.5 1 1 1];
W=[0 0 0 1 1 1];
% fold = '/Users/kiendljosef/Desktop/Forschung/solid_shell/Input_Files_edited/';
% name = strcat('SLo_2el_p2q2w2');
% CPm = textread(strcat(fold,name));
% CP = CP_matrix(p,q,r,U,V,W,CPm);
load CP_scord_comp_3CP.mat


[p,q,r,U,V,W,CP] = degree_elevate_solid(p,q,r,U,V,W,CP,deg-p,deg-q,deg-r);
R = refinement_vec(V,ref/2);
[CP,U,V,W] = knot_refine_solid(p,q,r,U,V,W,CP,R,R,[]);
nu = length(CP(:,1,1,1));
nv = length(CP(1,:,1,1));
nw = length(CP(1,1,:,1));

% material
E = 4.32e8;
nue = 0.0;
% Gauss Points
ngauss=[p+1 q+1 r+1];

% fixed supports
rb =[];
us=[0 0]; vs=[0 1]; ws=[0 1];  dirs=1;
rb = supports(rb,us,vs,ws,dirs,CP);
us=[0 1]; vs=[0 0]; ws=[0 1];  dirs=1;
rb = supports(rb,us,vs,ws,dirs,CP);
us=[0 1]; vs=[0 0]; ws=[0 1];  dirs=3;
rb = supports(rb,us,vs,ws,dirs,CP);
us=[0 1]; vs=[1 1]; ws=[0 1];  dirs=2;
rb = supports(rb,us,vs,ws,dirs,CP);

% load vector
f = -360;  fl=[];
ub=[0 1]; vb=[0 1]; wb=[0 1];  dirf=3;
fl = load_body(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,f,dirf);
% ub=[0 1]; vb=1; wb=[0 1];  dirf=2;  proj=0;
% fl = load_area(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,f,dirf,proj);
% ub=[0 1]; vb=0; wb=1;  dirf=3;  proj=0;
% fl = load_line(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,f,dirf,proj);

% compute displacement vector d and load vector including supports fs
d = solve(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);
%d = solve_ANS_standard(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);
%d = solve_ANS(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);
%d = solve_new(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);

CPd = CPresult(CP,d);

res(deg-1,iref)=(CPd(nu,nv,1,3)-CP(nu,nv,1,3)+CPd(nu,nv,nw,3)-CP(nu,nv,nw,3))/2;

end
end


% name2 = strcat('_deg',int2str(deg),'ref',int2str(ref));
% name3 = strcat(fold,name,name2,'.mat');
% save(name3,'en');

% visualize
% plot2in1(p,q,r,U,V,W,CP,CPd,rb,fl);
% plotNURBS_solid(p,q,r,U,V,W,CP);
% plotNURBS_solid(p,q,r,U,V,W,CPd);