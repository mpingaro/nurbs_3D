% 
% J. Kiendl
clear;

deg = 2;
degw = deg;
%refs = [2, 4, 8, 16];
refs = 8;
for iref = 1:length(refs)
ref = refs(iref);

clear CP;

% material
E = 4.32e8;
nue = 0.0;
% Gauss Points
ngauss=[3 3 3];

% geometry
R = 25;
L = 50;
t = 0.25;
alf = 40/180*pi;
x1=0; y1=0; z1=R; w1=1;
x2=R*tan(alf/2); y2=0; z2=R; w2=abs(cos(alf/2));
x3=R*sin(alf); y3=0; z3=R*cos(alf); w3=1;

p=2;
U=[0 0 0 1 1 1];
CP1 = [x1 y1 z1 w1; x2 y2 z2 w2; x3 y3 z3 w3];
% plotNURBS_full(U,CP1,p); axis equal;

% test circle
% Q = [R 0 0 1; R 0 R sqrt(0.5); 0 0 R 1];
% hold on; plotNURBS_full(U,Q,p); axis equal; hold off;

q=1; r=1;
V=[0 0 1 1]; W=[0 0 1 1];
CP(:,1,1,1) = CP1(:,1)/R*(R-t/2);
CP(:,1,1,2) = CP1(:,2);
CP(:,1,1,3) = CP1(:,3)/R*(R-t/2);
CP(:,1,1,4) = CP1(:,4);
CP(:,2,1,:) = CP(:,1,1,:);
CP(:,2,1,2) = CP(:,1,1,2)+L/2;
CP(:,:,2,:) = CP(:,:,1,:);
CP(:,:,2,1) = CP(:,:,1,1)/(R-t/2)*(R+t/2);
CP(:,:,2,3) = CP(:,:,1,3)/(R-t/2)*(R+t/2);

% plotNURBS_solid(p,q,r,U,V,W,CP)

[p,q,r,U,V,W,CP] = degree_elevate_solid(p,q,r,U,V,W,CP,deg-p,deg-q,degw-r);
Ru=refinement_vec(U,ref);
[CP,U,V,W] = knot_refine_solid(p,q,r,U,V,W,CP,Ru,Ru,[]);
% plotNURBS_solid(p,q,r,U,V,W,CP)

% diaphragm
rb =[];
us=[0 1]; vs=[0 0]; ws=[0 1];
rb = supports(rb,us,vs,ws,1,CP);
rb = supports(rb,us,vs,ws,3,CP);
% symmetry u=0
us=[0 0]; vs=[0 1]; ws=[0 1];
rb = supports(rb,us,vs,ws,1,CP);
% symmetry v=1
us=[0 1]; vs=[1 1]; ws=[0 1];
rb = supports(rb,us,vs,ws,2,CP);

% load vector
F = -360;  fl=[];
ub=[0 1]; vb=[0 1]; wb=[0 1];  dirf=3;
fl = load_body(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,F,dirf);
% fl2 = load_area([],ub,vb,1,p,q,r,U,V,W,CP,ngauss,-90,dirf,0);

% compute displacement vector d and load vector including supports fs
%d = solve(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);
%d = solve_ANS_standard(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);
d = solve_ANS(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);

CPd = CPresult(CP,d);
uhc = CPd(end,end,end,3)-CP(end,end,end,3);
resi(iref)= -uhc;

% visualize
% plot2in1(p,q,r,U,V,W,CP,CPda,rb,fl);
% plot1in1(p,q,r,U,V,W,CPdc,rb,fl);
% plot1in1(p,q,r,U,V,W,CPda,rb,fl);
% plot_strain_stress(p,q,r,U,V,W,CP,CPd,rb,fl,E,nue,d,2,10)
plot_strain_stress_1in1(p,q,r,U,V,W,CP,CPd,rb,fl,E,nue,d,2,6)

% uref = 0.3024;
% uhc = CPdc(end,end,end,3)-CP(end,end,end,3);
% uha = CPda(end,end,end,3)-CP(end,end,end,3);
% uhaxyz = CPdaxyz(end,end,end,3)-CP(end,end,end,3);
% resi=[-uhc -uha -uhaxyz uref]
% load 'res'; res(iref,:)=resi; save 'res' 'res';

end