function res = Example_study1_func(deg,degr,ref,t,lasd)

% NURBS parameters
p=1; q=1; r=1;
U=[0 0 1 1]; V=[0 0 1 1]; W=[0 0 1 1];
CP(:,:,1,1) = [0 0; 1 1];
CP(:,:,1,2) = [0 1; 0 1];
CP(:,:,1,3) = [0 0; 0 0];
CP(:,:,1,4) = [1 1; 1 1];
CP(:,:,2,:) = CP(:,:,1,:);
CP(:,:,2,3) = CP(:,:,1,3)+t;

% material
E = t^(-3);
nue = 0.2;

% REFINE
[p,q,r,U,V,W,CP] = degree_elevate_solid(p,q,r,U,V,W,CP,deg-p,deg-q,degr-r);
R = refinement_vec(U,ref);
[CP,U,V,W] = knot_refine_solid(p,q,r,U,V,W,CP,R,R,[]);
% plotNURBS_solid(p,q,r,U,V,W,CP)

% Gauss Points
ngauss=[p+1 q+1 r+1];

% fixed supports
rb =[];
us=[0 0]; vs=[0 1]; ws=[0 1];
rb = supports(rb,us,vs,ws,1,CP);
rb = supports(rb,us,vs,ws,2,CP);
rb = supports(rb,us,vs,ws,3,CP);
us=[1 1]; vs=[0 1]; ws=[0 1];
rb = supports(rb,us,vs,ws,1,CP);
rb = supports(rb,us,vs,ws,2,CP);
rb = supports(rb,us,vs,ws,3,CP);
us=[0 1]; vs=[0 0]; ws=[0 1]; 
rb = supports(rb,us,vs,ws,1,CP);
rb = supports(rb,us,vs,ws,2,CP);
rb = supports(rb,us,vs,ws,3,CP);
us=[0 1]; vs=[1 1]; ws=[0 1]; 
rb = supports(rb,us,vs,ws,1,CP);
rb = supports(rb,us,vs,ws,2,CP);
rb = supports(rb,us,vs,ws,3,CP);

% load vector
f = 1;  fl=[]; ub=[0 1]; vb=[0 1]; dirf=3;
if lasd == 1
  fl = load_area(fl,ub,vb,1,p,q,r,U,V,W,CP,ngauss,f,dirf,0);
elseif lasd == 2
  fl = load_body(fl,ub,vb,[0 1],p,q,r,U,V,W,CP,ngauss,f/t,dirf);
end

% compute displacement vector d and load vector including supports fs
d = solve2(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);
% d = solve_ANS2(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);
CPd = CPresult(CP,d);

res.U=U;
res.p=p;
res.r=r;
res.d=d;
res.CP=CP;
res.CPd=CPd;

