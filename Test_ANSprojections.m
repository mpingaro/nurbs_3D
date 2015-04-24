clear;

a = [-sqrt(1/3) sqrt(1/3)];
ua = a/2+.5;
b = [-sqrt(3/5) 0 sqrt(3/5)];
ub = b/2+.5;

% Test for M1 =========================================
p=1;q=2;U=[0 0 1 1];V=[0 0 0 1 1 1];
CP(:,:,1)=[-1 -1 -1; 1 1 1];
CP(:,:,2)=[-1 0 1; -1 0 1];
CP(:,:,3)=[-0.5 0 0; 0 0 0.5];
CP(:,:,4)=[1 1 1; 1 1 1];
plotNURBS_surf(p,q,U,V,CP); hold on;
create_el_edges(p,q,U,V,CP); axis equal; %hold off;
k = 1;
for j=1:3
  for i=1:2
    Q(k,:) = get_point_surf(p,0,ua(i),U,q,0,ub(j),V,CP);
    k = k+1;
  end
end
ik = 1;
for j=1:3
  for i=1:3
    P(ik,:) = get_point_surf(p,0,ub(i),U,q,0,ub(j),V,CP);
    ik = ik+1;
  end
end
for k=1:6
  plot3(Q(k,1),Q(k,2),Q(k,3),'^');
end
for k=1:9
  plot3(P(k,1),P(k,2),P(k,3),'o');
end
M1 = [1/2+sqrt(9/20)  1/2-sqrt(9/20) 0 0 0 0;...
      1/2             1/2            0 0 0 0;...
      1/2-sqrt(9/20)  1/2+sqrt(9/20) 0 0 0 0;...
      0 0 1/2+sqrt(9/20)  1/2-sqrt(9/20) 0 0;...
      0 0 1/2             1/2            0 0;...
      0 0 1/2-sqrt(9/20)  1/2+sqrt(9/20) 0 0;...
      0 0 0 0 1/2+sqrt(9/20)  1/2-sqrt(9/20);...
      0 0 0 0 1/2             1/2           ;...
      0 0 0 0 1/2-sqrt(9/20)  1/2+sqrt(9/20)];
        
P2 = M1*Q;
for k=1:9
  plot3(P2(k,1),P2(k,2),P2(k,3),'x');
end
% P2-P
axis equal;
view(2);
hold off;
% Test for M2 =========================================
clear 'CP' 'P' 'Q' 'P2';
p=2;q=1;U=[0 0 0 1 1 1];V=[0 0 1 1];
CP(:,:,1)=[-1 -1; 0 0; 1 1];
CP(:,:,2)=[-1 1; -1 1; -1 1];
CP(:,:,3)=[-0.5 0; 0 0; 0 0.5];
CP(:,:,4)=[1 1; 1 1; 1 1];
plotNURBS_surf(p,q,U,V,CP); hold on;
create_el_edges(p,q,U,V,CP); axis equal; %hold off;
k = 1;
for j=1:2
  for i=1:3
    Q(k,:) = get_point_surf(p,0,ub(i),U,q,0,ua(j),V,CP);
    k = k+1;
  end
end
ik = 1;
for j=1:3
  for i=1:3
    P(ik,:) = get_point_surf(p,0,ub(i),U,q,0,ub(j),V,CP);
    ik = ik+1;
  end
end
for k=1:6
  plot3(Q(k,1),Q(k,2),Q(k,3),'^');
end
for k=1:9
  plot3(P(k,1),P(k,2),P(k,3),'o');
end
M2 = [1/2+sqrt(9/20) 0 0 1/2-sqrt(9/20) 0 0;...
      0 1/2+sqrt(9/20) 0 0 1/2-sqrt(9/20) 0;...
      0 0 1/2+sqrt(9/20) 0 0 1/2-sqrt(9/20);...
      1/2            0 0 1/2            0 0;...
      0 1/2            0 0 1/2            0;...
      0 0 1/2            0 0 1/2           ;...
      1/2-sqrt(9/20) 0 0 1/2+sqrt(9/20) 0 0;...
      0 1/2-sqrt(9/20) 0 0 1/2+sqrt(9/20) 0;...
      0 0 1/2-sqrt(9/20) 0 0 1/2+sqrt(9/20)];
        
P2 = M2*Q;
for k=1:9
  plot3(P2(k,1),P2(k,2),P2(k,3),'x');
end
% P2-P
axis equal;
view(2);
hold off;
% Test for M3 =========================================
clear 'CP' 'P' 'Q' 'P2';
p=1;q=1;U=[0 0 1 1];V=[0 0 1 1];
CP(:,:,1)=[-1 -1; 1 1];
CP(:,:,2)=[-1 1; -1 1];
CP(:,:,3)=[-0.2 0.3; 0.5 -.3];
CP(:,:,4)=[1 1; 1 1];
plotNURBS_surf(p,q,U,V,CP); hold on;
create_el_edges(p,q,U,V,CP); axis equal; %hold off;
k = 1;
for j=1:2
  for i=1:2
    Q(k,:) = get_point_surf(p,0,ua(i),U,q,0,ua(j),V,CP);
    k = k+1;
  end
end
ik = 1;
for j=1:3
  for i=1:3
    P(ik,:) = get_point_surf(p,0,ub(i),U,q,0,ub(j),V,CP);
    ik = ik+1;
  end
end
for k=1:4
  plot3(Q(k,1),Q(k,2),Q(k,3),'^');
end
for k=1:9
  plot3(P(k,1),P(k,2),P(k,3),'o');
end
M3 = 1/4*[14/5+sqrt(36/5) -4/5          -4/5          14/5-sqrt(36/5);...
          1+sqrt(9/5)   1+sqrt(9/5)   1-sqrt(9/5)   1-sqrt(9/5);...
          -4/5          14/5+sqrt(36/5) 14/5-sqrt(36/5) -4/5;...
          1+sqrt(9/5)   1-sqrt(9/5)   1+sqrt(9/5)   1-sqrt(9/5);...
          1           1           1           1;...
          1-sqrt(9/5)   1+sqrt(9/5)   1-sqrt(9/5)   1+sqrt(9/5);...
          -4/5          14/5-sqrt(36/5) 14/5+sqrt(36/5) -4/5;...
          1-sqrt(9/5)   1-sqrt(9/5)   1+sqrt(9/5)   1+sqrt(9/5);...
          14/5-sqrt(36/5) -4/5          -4/5          14/5+sqrt(36/5)];
        
P2 = M3*Q;
for k=1:9
  plot3(P2(k,1),P2(k,2),P2(k,3),'x');
end
% P2-P
axis equal;
view(2);
hold off;