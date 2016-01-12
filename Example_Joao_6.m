% Example from ANS solid-shell paper
% 6. Pinched Cylinder
clear;
clc;

deg = 3;
% refs = [2, 4, 8, 16, 32];
refs = 16;

% res = zeros(4,5,2);

%for deg = 3:3
for iref = 1:length(refs)
  ref = refs(iref);

clear CP;
  
p=2;
q=1;
r=1;
U=[0 0 0 1 1 1];
V=[0 0 1 1];
W=[0 0 1 1];
w=sqrt(0.5);
t = 3;
% Control Point coordinates
CP(:,:,1,1)=[1 1; 1 1; 0 0]*(300-t/2);
CP(:,:,1,2)=[0 0; 1 1; 1 1]*(300-t/2);
CP(:,:,1,3)=[0 1; 0 1; 0 1]*300;
CP(:,:,1,4)=[1 1; w w; 1 1];
CP(:,:,2,1)=[1 1; 1 1; 0 0]*(300+t/2);
CP(:,:,2,2)=[0 0; 1 1; 1 1]*(300+t/2);
CP(:,:,2,3:4)=CP(:,:,1,3:4);

[p,q,r,U,V,W,CP] = degree_elevate_solid(p,q,r,U,V,W,CP,deg-p,deg-q,deg-r);
R = refinement_vec(V,ref);
[CP,U,V,W] = knot_refine_solid(p,q,r,U,V,W,CP,R,R,[]);
nu = length(CP(:,1,1,1));
nv = length(CP(1,:,1,1));
nw = length(CP(1,1,:,1));

% material
E = 3e6;
nue = 0.3;
% Gauss Points
ngauss=[p+1 q+1 r+1];

% fixed supports
rb =[];
us=[0 1]; vs=[0 0]; ws=[0 1];  dirs=1;
rb = supports(rb,us,vs,ws,dirs,CP);
us=[0 1]; vs=[0 0]; ws=[0 1];  dirs=2;
rb = supports(rb,us,vs,ws,dirs,CP);
us=[0 1]; vs=[1 1]; ws=[0 1];  dirs=3;
rb = supports(rb,us,vs,ws,dirs,CP);
us=[0 0]; vs=[0 1]; ws=[0 1];  dirs=2;
rb = supports(rb,us,vs,ws,dirs,CP);
us=[1 1]; vs=[0 1]; ws=[0 1];  dirs=1;
rb = supports(rb,us,vs,ws,dirs,CP);

% load vector
f = -0.25/t;  fl=[];
ub=0; vb=1; wb=[0 1];  dirf=1;  proj=0;
fl = load_line(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,f,dirf,proj);

% compute displacement vector d and load vector including supports fs
d = solve(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);
%d = solve_ANS_standard(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);
%d = solve_ANS(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);

CPd = CPresult(CP,d);

% d2=200*d;
% CPd2 = CPresult(CP,d2);
% plot2in1(p,q,r,U,V,W,CP,CPd2,rb,fl);

S = get_point_solid(0,1,0,0,0,0,p,q,r,U,V,W,CP);
Sd = get_point_solid(0,1,0,0,0,0,p,q,r,U,V,W,CPd);
disp1x = Sd(1)-S(1);
S = get_point_solid(0,1,1,0,0,0,p,q,r,U,V,W,CP);
Sd = get_point_solid(0,1,1,0,0,0,p,q,r,U,V,W,CPd);
disp2x = Sd(1)-S(1);
res(deg-1,iref)=(disp1x+disp2x)/2;

% res(1:nw,iref,deg-1)=CPd(1,1,1:nw,1)-CP(1,1,1:nw,1);

% visualize
plot_strain_stress_1in1(p,q,r,U,V,W,CP,CPd,rb,fl,E,nue,d,2,6)

end
%end


% name2 = strcat('_deg',int2str(deg),'ref',int2str(ref));
% name3 = strcat(fold,name,name2,'.mat');
% save(name3,'en');

% visualize
% plot2in1(p,q,r,U,V,W,CP,CPd,rb,fl);
% plotNURBS_solid(p,q,r,U,V,W,CP);
% plotNURBS_solid(p,q,r,U,V,W,CPd);