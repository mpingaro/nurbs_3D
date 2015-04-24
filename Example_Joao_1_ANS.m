% Example from ANS solid-shell paper
% 1+2. cantilever beam with mesh refinement or thickness change
% with p=2 and ANS
clear;

deg = 2;
%refs = [2 4 8 16];
%ts = [5 1 0.5 0.1 0.05 0.02 0.01];
refs = [8];
ts = [5];

for it = 1:length(ts)-2
  t=ts(it);
for iref=3:3%1:length(refs)
  ref=refs(iref);
  
p=1;
q=1;
r=1;
U=[0 0 1 1];
V=[0 0 1 1];
W=[0 0 1 1];
clear CP;
% Control Point coordinates
CP(:,:,1,1)=[0 0; 1 1]*100;
CP(:,:,1,2)=[0 1; 0 1];
CP(:,:,1,3)=[0 0; 0 0];
CP(:,:,1,4)=[1 1; 1 1];
CP(:,:,2,:)=CP(:,:,1,:);
CP(:,:,2,3)=CP(:,:,1,3)+t;
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
f = 0.1*t^(1/2);  fl=[];
% ub=[0 1]; vb=[0 1]; wb=[0 1];  dirf=3;
% fl = load_body(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,f,dirf);
ub=1; vb=[0 1]; wb=[0 1];  dirf=3;  proj=0;
fl = load_area(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,f,dirf,proj);
% ub=[0 1]; vb=0; wb=1;  dirf=3;  proj=0;
% fl = load_line(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,f,dirf,proj);

% compute displacement vector d and load vector including supports fs
[dc,fs] = solve(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);
[da,fs] = solve_ANS(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);
CPdc = CPresult(CP,dc);
CPda = CPresult(CP,da);

enc = 0.5*fl'*dc;
ena = 0.5*fl'*da;

% res(iref) = en;
resc(it,:) = [t,enc/20];
resa(it,:) = [t,ena/20];

end
end

sl=100./resc(:,1);
semilogx(sl,resc(:,2),'-x'); 
hold on;
semilogx(sl,resa(:,2),'-o');
ylim([0.6 1.1]);
grid on;
legend('p=2','p=2 ANS xyz');
hold off;

% name2 = strcat('_deg',int2str(deg),'ref',int2str(ref));
% name3 = strcat(fold,name,name2,'.mat');
% save(name3,'en');

% visualize
plot2in1(p,q,r,U,V,W,CP,CPda,rb,fl)