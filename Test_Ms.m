% compare my transformation matrices Ms with Robertts N*M^-1

clear;
load Ms; % my matrices
p1=1; U1=[0 0 1 1];     % linear 
p2=2; U2=[0 0 0 1 1 1]; % quadratic
% Control points are not really used. But for the Nurbs basis functions routine, 
% you need to provide them for the control weights. So, I just set everything =1.
CPn =ones(3,3,4); % 3x3 integration pts
CPm1=ones(2,3,4); % 2x3 tying pts for xixi
CPm2=ones(3,2,4); % 3x2 tying pts for etaeta
CPm3=ones(2,2,4); % 2x2 tying pts for xieta
u1=[-sqrt(1/3) sqrt(1/3)]/2+0.5;    % linear Gauss coordinates
u2=[-sqrt(3/5) 0 sqrt(3/5)]/2+0.5;  % quadratic Gauss coordinates
% M_xixi ---------------------
k=1; 
for j=1:3 
  for i=1:3 
    [R,dR] = deriv_Nurbsbasisfunc01(p1,0,u2(i),U1,p2,0,u2(j),U2,CPn); 
    N1(k,:)=R; k=k+1; 
  end 
end
k=1; 
for j=1:3 
  for i=1:2 
    [R,dR] = deriv_Nurbsbasisfunc01(p1,0,u1(i),U1,p2,0,u2(j),U2,CPm1); 
    M1r(k,:)=R; k=k+1; 
  end 
end
M1p = N1*inv(M1r);
dM1=M1p-M1 
% M_etaeta ---------------------
k=1; 
for j=1:3 
  for i=1:3 
    [R,dR] = deriv_Nurbsbasisfunc01(p2,0,u2(i),U2,p1,0,u2(j),U1,CPn); 
    N2(k,:)=R; k=k+1; 
  end 
end
k=1; 
for j=1:2 
  for i=1:3 
    [R,dR] = deriv_Nurbsbasisfunc01(p2,0,u2(i),U2,p1,0,u1(j),U1,CPm2); 
    M2r(k,:)=R; k=k+1; 
  end 
end
M2p = N2*inv(M2r);
dM2=M2p-M2 
% M_xieta ---------------------
k=1; 
for j=1:3 
  for i=1:3 
    [R,dR] = deriv_Nurbsbasisfunc01(p1,0,u2(i),U1,p1,0,u2(j),U1,CPn); 
    N3(k,:)=R; k=k+1; 
  end 
end
k=1; 
for j=1:2 
  for i=1:2 
    [R,dR] = deriv_Nurbsbasisfunc01(p1,0,u1(i),U1,p1,0,u1(j),U1,CPm3); 
    M3r(k,:)=R; k=k+1; 
  end 
end
M3p = N3*inv(M3r);
dM3=M3p-M3 