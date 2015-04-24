function plot1in1(p,q,r,U,V,W,CP,rb,f)
% including boundary conditions and loads

% calculate plot points for:
% 1. original geometry, supports and load arrows
[xs1,ys1,zs1] = create_supports(CP,rb);
[xf1,yf1,zf1] = create_arrows(CP,f);


% geometry
plotNURBS_solid(p,q,r,U,V,W,CP);
hold on;

% supports
for k =1:length(xs1(:,1))
  plot3(xs1(k,:),ys1(k,:),zs1(k,:),'Linewidth',2,'Color','black');
end

% load arrows
for k =1:length(xf1(:,1))
  plot3(xf1(k,:),yf1(k,:),zf1(k,:),'Linewidth',5);
  plot3(xf1(k,1),yf1(k,1),zf1(k,1),'Marker','d','MarkerFaceColor','blue','MarkerSize',10);
end

axis equal;
camlight left;
hold off;