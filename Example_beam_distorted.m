% Example from ANS solid-shell paper
% cantilever beam with mesh distorted refinement and fix thickness
clear all; close all; clc;

deg = 2;
%refs = [2 4 8 16];
refs = 8;

sl = [10, 20, 100, 200, 1000, 2000, 5000, 10000];
w = 1;
L = 100;
res = zeros(length(refs),2);
ang = pi/6;
ds = 1/2*tan(ang);

for isl=1:length(sl)
t=L/sl(isl);
for iref=1:length(refs)
ref=refs(iref);

p=2;
q=2;
r=2;
U=[0 0 0 1 1 1];
V=[0 0 0 1 1 1];
W=[0 0 0 1 1 1];
clear CP;
% % Control points and weight
coefs(:,:,1,1) = [0, L/2+ds, L;
                  0, 0, 0;
                  0, 0, 0;
                  1, 1, 1];
coefs(:,:,2,1) = [0, L/2, L;
                  w/2, w/2, w/2;
                  0, 0, 0;
                  1, 1, 1];
coefs(:,:,3,1) = [0, L/2-ds, L;
                  w, w, w;
                  0, 0, 0;
                  1, 1, 1];                  

coefs(:,:,1,2) = [0, L/2+ds, L;
                  0, 0, 0;
                  t/2, t/2, t/2;
                  1, 1, 1];
coefs(:,:,2,2) = [0, L/2, L;
                  w/2, w/2, w/2;
                  t/2, t/2, t/2;
                  1, 1, 1];
coefs(:,:,3,2) = [0, L/2-ds, L;
                  w, w, w;
                  t/2, t/2, t/2;
                  1, 1, 1]; 

coefs(:,:,1,3) = [0, L/2+ds, L;
                  0, 0, 0;
                  t, t, t;
                  1, 1, 1];
coefs(:,:,2,3) = [0, L/2, L;
                  w/2, w/2, w/2;
                  t, t, t;
                  1, 1, 1];
coefs(:,:,3,3) = [0, L/2-ds, L;
                  w, w, w;
                  t, t, t;
                  1, 1, 1]; 
              
for i=1:3
    for j=1:3
        for h=1:3
            for k=1:4
                CP(i,j,h,k) = coefs(k,i,j,h);
            end
        end
    end
end
% % Control Point coordinates
% CP(:,:,1,1)=[0 0; 1 1]*L;
% CP(:,:,1,2)=[0 1; 0 1];
% CP(:,:,1,3)=[0 0; 0 0];
% CP(:,:,1,4)=[1 1; 1 1];
% CP(:,:,2,:)=CP(:,:,1,:);
% CP(:,:,2,3)=CP(:,:,1,3)+t;
% [p,q,r,U,V,W,CP] = degree_elevate_solid(p,q,r,U,V,W,CP,deg-p,deg-q,deg-r);
%plotNURBS_solid(p,q,r,U,V,W,CP);

R = refinement_vec(U,ref);
[CP,U,V,W] = knot_refine_solid(p,q,r,U,V,W,CP,R,[],[]);
nu=length(CP(:,1,1,1)); nv=length(CP(1,:,1,1)); nw=length(CP(1,1,:,1));
% Plot geometry
% figure, plotNURBS_solid(p,q,r,U,V,W,CP);

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
f = E/4/t*(t/L)^3;  fl=[];
% ub=[0 1]; vb=[0 1]; wb=[0 1];  dirf=3;
% fl = load_body(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,f,dirf);
ub=1; vb=[0 1]; wb=[0 1];  dirf=3;  proj=0;
fl = load_area(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,f,dirf,proj);
% ub=[0 1]; vb=0; wb=1;  dirf=3;  proj=0;
% fl = load_line(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,f,dirf,proj);

% compute displacement vector d and load vector including supports fs
[d,fs] = solve_ANS_standard(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);
CPd = CPresult(CP,d);

S = get_point_solid(1,1,1,0,0,0,p,q,r,U,V,W,CP);
Sd = get_point_solid(1,1,1,0,0,0,p,q,r,U,V,W,CPd);
dispz = Sd(3)-S(3);

%res(iref,:) = [ref,dispz];
res(isl,:) = [sl(isl), dispz];
end
end
