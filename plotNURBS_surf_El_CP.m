function plotNURBS_surf_El_CP(p,q,U,V,CP)
% plots the surface, elements and control points

[X,Y,Z] = create_surf(p,q,U,V,CP);

% geometry
surf(X,Y,Z,'FaceColor','green','EdgeColor','none');
hold on;

% element edges
create_el_edges(p,q,U,V,CP)

% control points and polygon
create_conpolygon(CP)

% camlight left; lighting phong;
axis equal;
xlabel('');
ylabel('');
zlabel('');
hold off;
% view(2);