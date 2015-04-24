% corrected vertical cantilever
% J. Kiendl
clear;

refs = 2;
% refs = [2 4 8 16 32];
L=100; t=0.1; w=10;

% material
% t = 1;
E = 1e6;
nue = 0.0;

% Gauss Points
ngauss=[3 3 3];


for iref=1:length(refs)

  p=1; q=1; r=1;
  U=[0 0 1 1]; V=[0 0 1 1]; W=[0 0 1 1];
  clear CP;
  CP(:,:,1,1) = [0 0; 0 0];
  CP(:,:,1,2) = [0 1; 0 1]*w;
  CP(:,:,1,3) = [0 0; 1 1]*L;
  CP(:,:,1,4) = [1 1; 1 1];
  CP(:,:,2,:) = CP(:,:,1,:);
  CP(:,:,2,1) = CP(:,:,1,1)-t;
%   plotNURBS_solid(p,q,r,U,V,W,CP)

  [p,q,r,U,V,W,CP] = degree_elevate_solid(p,q,r,U,V,W,CP,1,1,1);
  Ru=refinement_vec(U,refs(iref));
  [CP,U,V,W] = knot_refine_solid(p,q,r,U,V,W,CP,Ru,[0.5],[]);
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
  F = -1;  fl=[];
  ub=1; vb=[0 1]; wb=[0 1];  dirf=1;
  fl = load_area(fl,ub,vb,wb,p,q,r,U,V,W,CP,ngauss,F,dirf,0);
  % fl2 = load_area(fl2,ub,vb,wb,p2,q2,r2,U2,V2,W2,CP2,ngauss,-fq,dirf,proj);
  
  % compute displacement vector d and load vector including supports fs
  [dc,fs] = solve(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);
  [da,fs] = solve_ANS(p,U,q,V,r,W,CP,E,nue,ngauss,fl,rb);
  CPdc = CPresult(CP,dc);
  CPda = CPresult(CP,da);
  
  % visualize
%   plot2in1(p,q,r,U,V,W,CP,CPda,rb,fl);
  % plot_strain_stress(p,q,r,U,V,W,CP,CPd,rb,fl,E,nue,d,2,2)
  
  enhc = 0.5*dc'*fl;
  enha(iref) = 0.5*da'*fl;
%   enh2 = 0.5*d'*fs
  ene = abs(2*(F*w*t)^2*L^3/E/w/t^3);
%   eneT = abs(2*(F*w*t)^2*L^3/E/w/t^3 + F^2*6*(1+nue)*L/(5*E*w*t));
%   err(iref)=(enha(iref)-ene)/ene
%   errT(iref)=(enha(iref)-eneT)/eneT
  [enhc enha ene]/ene

end