% Example from ANS solid-shell paper
% 3. circular beam with thickness change
clear;

deg = 3;
ref = 10;
ts = [1 0.5 0.1 0.05 0.01 0.005 0.001];

for it = 1:length(ts)
  t=ts(it);
  
p=2;
q=1;
r=1;
U=[0 0 0 1 1 1];
V=[0 0 1 1];
W=[0 0 1 1];
clear CP;
w=sqrt(1/2);
% Control Point coordinates
CP(:,1,1,1)=[0; 1; 1]*(10-t/2);
CP(:,1,1,2)=[0; 0; 0];
CP(:,1,1,3)=[1; 1; 0]*(10-t/2);
CP(:,1,1,4)=[1; w; 1];
CP(:,1,2,1:3)=CP(:,1,1,1:3)/(10-t/2)*(10+t/2);
CP(:,1,2,4)=CP(:,1,1,4);
CP(:,2,:,:)=CP(:,1,:,:);
CP(:,2,:,2)=CP(:,1,:,2)+1;
% plotNURBS_solid(p,q,r,U,V,W,CP);

[p,q,r,U,V,W,CP] = degree_elevate_solid(p,q,r,U,V,W,CP,deg-p,deg-q,deg-r);
R = refinement_vec(U,ref);
[CP,U,V,W] = knot_refine_solid(p,q,r,U,V,W,CP,R,[],[]);
nu=length(CP(:,1,1,1)); nv=length(CP(1,:,1,1)); nw=length(CP(1,1,:,1));
% plotNURBS_solid(p,q,r,U,V,W,CP);

% material
E = 1e3;
nue = 0.0;
% Gauss Points
ngauss=[p+1 q+1 r+1];

% fixed supports
rb =[];
us=[0 0]; vs=[0 1]; ws=[0 1];  dirs=1;
rb = supports(rb,us,vs,ws,dirs,CP);
us=[0 0]; vs=[0 1]; ws=[0 1];  dirs=2;
rb = supports(rb,us,vs,ws,dirs,CP);
us=[0 0]; vs=[0 1]; ws=[0 1];  dirs=3;
rb = supports(rb,us,vs,ws,dirs,CP);

% load vector
f = -0.1*t^2;  fl=[];
% ub=[0 1]; vb=[0 1]; wb=[0 1];  dirf=3;
% fl = load_body(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,f,dirf);
ub=1; vb=[0 1]; wb=[0 1];  dirf=1;  proj=0;
fl = load_area(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,f,dirf,proj);
% ub=[0 1]; vb=0; wb=1;  dirf=3;  proj=0;
% fl = load_line(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,f,dirf,proj);

% compute displacement vector d and load vector including supports fs
[d,fs] = solve(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);
CPd = CPresult(CP,d);

S = get_point_solid(1,0,0,0,0,0,p,q,r,U,V,W,CP);
Sd = get_point_solid(1,0,0,0,0,0,p,q,r,U,V,W,CPd);
disp1x = Sd(1)-S(1);
S = get_point_solid(1,0,1,0,0,0,p,q,r,U,V,W,CP);
Sd = get_point_solid(1,0,1,0,0,0,p,q,r,U,V,W,CPd);
disp2x = Sd(1)-S(1);
dispx=(disp1x+disp2x)/2;

res(it,:) = [t,dispx];

end

res

% name2 = strcat('_deg',int2str(deg),'ref',int2str(ref));
% name3 = strcat(fold,name,name2,'.mat');
% save(name3,'en');

% visualize
% plot2in1(p,q,r,U,V,W,CP,CPd,rb,fl)