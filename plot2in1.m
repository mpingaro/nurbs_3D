function plot2in1(p,q,r,U,V,W,CP1,CP2,rb,f)
% plots the original and the deformed geometry in one figure,
% including boundary conditions and loads

% calculate plot points for:
% 1. original geometry, supports and load arrows
[xs1,ys1,zs1] = create_supports(CP1,rb);
[xf1,yf1,zf1] = create_arrows(CP1,f);
% 2. deformed geometry, supports and load arrows
[xs2,ys2,zs2] = create_supports(CP2,rb);
[xf2,yf2,zf2] = create_arrows(CP2,f);

% run plot two times. first run only for getting the plot limits
for s = 1:2
  
  % FIRST WINDOW: UNDEFORMED
  subplot(2,1,1);
  
  % geometry
  plotNURBS_solid(p,q,r,U,V,W,CP1);
  hold;
  
  % supports
  for k =1:length(xs1(:,1))
    plot3(xs1(k,:),ys1(k,:),zs1(k,:),'Linewidth',2,'Color','black');
  end
  
  % load arrows
  for k =1:length(xf1(:,1))
    plot3(xf1(k,:),yf1(k,:),zf1(k,:),'Linewidth',5);
    plot3(xf1(k,1),yf1(k,1),zf1(k,1),'Marker','d','MarkerFaceColor','blue','MarkerSize',10);
  end

  % scaling plot
  if (s==1)
    axis equal;
    v1 = axis;
  elseif (s==2)
    axis equal;
    camlight left;
    xlim([x1 x2]);
    ylim([y1 y2]);
    zlim([z1 z2]);
    xlabel('x','FontSize',14);
    ylabel('y','FontSize',14);
    zlabel('z','FontSize',14);
    title ('original geometry');
  end
  hold off;

  
  % SECOND WINDOW: DEFORMED
  subplot(2,1,2);
  
  % geometry
  plotNURBS_solid(p,q,r,U,V,W,CP2);
  hold;
  
  %supports
  for k =1:length(xs2(:,1))
    plot3(xs2(k,:),ys2(k,:),zs2(k,:),'Linewidth',2,'Color','black');
  end
  
  % load arrows  
  for k =1:length(xf2(:,1))
    plot3(xf2(k,:),yf2(k,:),zf2(k,:),'Linewidth',5);
    plot3(xf2(k,1),yf2(k,1),zf2(k,1),'Marker','d','MarkerFaceColor','blue','MarkerSize',10);
  end
  
  % scaling plot
  if (s==1)
    axis equal;
    v2 = axis;
    x1 = min(v1(1),v2(1));
    x2 = max(v1(2),v2(2));
    y1 = min(v1(3),v2(3));
    y2 = max(v1(4),v2(4));
    z1 = min(v1(5),v2(5));
    z2 = max(v1(6),v2(6));
  elseif (s==2)
    axis equal;
    camlight left;
    xlim([x1 x2]);
    ylim([y1 y2]);
    zlim([z1 z2]);
    xlabel('x','FontSize',14);
    ylabel('y','FontSize',14);
    zlabel('z','FontSize',14);
    title ('deformed geometry');
  end  
  hold off;
end