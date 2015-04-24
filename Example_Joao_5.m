% Example from ANS solid-shell paper
% 5. Hemisphere
clear;

% deg = 2;
refs = [2 4 8 16];
% res = zeros(4,5,2);

for deg = 3:3
for iref = 1:length(refs)
  ref = refs(iref);

clear CP;
  
p=2;
q=2;
r=2;
U=[0 0 0 1 1 1];
V=[0 0 0 1 1 1];
W=[0 0 0 1 1 1];
w=sqrt(0.5);
% Control Point coordinates
CP(:,:,2,1)=[1 1 0; 1 1 0; 0 0 0]*10;
CP(:,:,2,2)=[0 0 0; 1 1 0; 1 1 0]*10;
CP(:,:,2,3)=[0 1 1; 0 1 1; 0 1 1]*10;
CP(:,:,2,4)=[1 w 1; w 0.5 w; 1 w 1];
CP(:,:,1,1:3)=CP(:,:,2,1:3)*0.998; CP(:,:,1,4)=CP(:,:,2,4);
CP(:,:,3,1:3)=CP(:,:,2,1:3)*1.002; CP(:,:,3,4)=CP(:,:,2,4);

[p,q,r,U,V,W,CP] = degree_elevate_solid(p,q,r,U,V,W,CP,deg-p,deg-q,deg-r);
R = refinement_vec(V,ref);
[CP,U,V,W] = knot_refine_solid(p,q,r,U,V,W,CP,R,R,[]);
nu = length(CP(:,1,1,1));
nv = length(CP(1,:,1,1));
nw = length(CP(1,1,:,1));

% material
E = 6.825e7;
nue = 0.3;
% Gauss Points
ngauss=[p+1 q+1 r+1];

% fixed supports
rb =[];
us=[0 0]; vs=[0 1]; ws=[0 1];  dirs=2;
rb = supports(rb,us,vs,ws,dirs,CP);
us=[1 1]; vs=[0 1]; ws=[0 1];  dirs=1;
rb = supports(rb,us,vs,ws,dirs,CP);
us=[0 1]; vs=[1 1]; ws=[0 1];  dirs=1;
rb = supports(rb,us,vs,ws,dirs,CP);
us=[0 1]; vs=[1 1]; ws=[0 1];  dirs=2;
rb = supports(rb,us,vs,ws,dirs,CP);
us=[0 1]; vs=[1 1]; ws=[0 1];  dirs=3;
rb = supports(rb,us,vs,ws,dirs,CP);

% load vector
f = 1/0.04;  fl=[];
% ub=[0 1]; vb=[0 1]; wb=[0 1];  dirf=3;
% fl = load_body(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,f,dirf);
% ub=[0 1]; vb=1; wb=[0 1];  dirf=2;  proj=0;
% fl = load_area(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,f,dirf,proj);
ub=0; vb=0; wb=[0 1];  dirf=1;  proj=0;
fl = load_line(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,f,dirf,proj);
ub=1; vb=0; wb=[0 1];  dirf=2;  proj=0;
fl = load_line(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,-f,dirf,proj);

% compute displacement vector d and load vector including supports fs
d = solve_new(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);
CPd = CPresult(CP,d);

% d2=200*d;
% CPd2 = CPresult(CP,d2);
% plot2in1(p,q,r,U,V,W,CP,CPd2,rb,fl);

S = get_point_solid(0,0,0,0,0,0,p,q,r,U,V,W,CP);
Sd = get_point_solid(0,0,0,0,0,0,p,q,r,U,V,W,CPd);
disp1x = Sd(1)-S(1);
S = get_point_solid(0,0,1,0,0,0,p,q,r,U,V,W,CP);
Sd = get_point_solid(0,0,1,0,0,0,p,q,r,U,V,W,CPd);
disp2x = Sd(1)-S(1);
res(deg-1,iref)=(disp1x+disp2x)/2;

% res(1:nw,iref,deg-1)=CPd(1,1,1:nw,1)-CP(1,1,1:nw,1);

end
end


% name2 = strcat('_deg',int2str(deg),'ref',int2str(ref));
% name3 = strcat(fold,name,name2,'.mat');
% save(name3,'en');

% visualize
% plot2in1(p,q,r,U,V,W,CP,CPd,rb,fl);
% plotNURBS_solid(p,q,r,U,V,W,CP);
% plotNURBS_solid(p,q,r,U,V,W,CPd);