function plotNURBS_full(U,CP,p)
% plots curve, control points and control polygon in x-z-plane

m=length(U);
n=length(CP(:,1));
% check_input(p,m,n,0,1,0)

eps=10e-10;

% CALCULATE POINTS
% curve
k=1;                %counting index for points
s=(U(m)-U(1))/1000;  %incremental step for u
u=U(1);
while u <= U(m)+eps
  i = findspan(u,U,n);   
  Cu(k,1:3) = get_point_curv(p,i,u,U,CP);  
  k=k+1;
  u=u+s;
end

% knots
for i2 = p+2:m-p-1
  i = findspan(U(i2),U,n);
  k=i2-p-1; 
  Knots(k,1:3) = get_point_curv(p,i,U(i2),U,CP);
end

% PLOT
% curve
plot (Cu(:,1),Cu(:,3),'r');
xlabel('x','FontSize',14);
ylabel('z','FontSize',14);
axis equal;
hold on;
% Knots
if (m>2*(p+1))
  plot (Knots(:,1),Knots(:,3),'linestyle','none','marker','x','color','b','markersize',12);
end
% Control Points
% plot (CP(:,1),CP(:,3),'marker','o','color','red');
plot (CP(:,1),CP(:,3),'--ok');
hold off;

