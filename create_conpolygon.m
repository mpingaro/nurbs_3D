function create_conpolygon(CP)
% plots control points and polygon

for k = 1:length(CP(1,1,:,1))-1
  for j = 1:length(CP(1,:,1,1))-1
    for i = 1:length(CP(:,1,1,1))-1
      plot3(CP(i:i+1,j,k,1),CP(i:i+1,j,k,2),CP(i:i+1,j,k,3),'--or');
      plot3(CP(i,j:j+1,k,1),CP(i,j:j+1,k,2),CP(i,j:j+1,k,3),'--or');
      CPz(:,1:3)=CP(i,j,:,1:3);
      plot3(CPz(k:k+1,1),CPz(k:k+1,2),CPz(k:k+1,3),'--or');
    end
    i=i+1;
    plot3(CP(i,j:j+1,k,1),CP(i,j:j+1,k,2),CP(i,j:j+1,k,3),'--or');
    CPz(:,1:3)=CP(i,j,:,1:3);
    plot3(CPz(k:k+1,1),CPz(k:k+1,2),CPz(k:k+1,3),'--or');
  end
  j=j+1;
  for i = 1:length(CP(:,1,1,1))-1
    plot3(CP(i:i+1,j,k,1),CP(i:i+1,j,k,2),CP(i:i+1,j,k,3),'--or');
    CPz(:,1:3)=CP(i,j,:,1:3);
    plot3(CPz(k:k+1,1),CPz(k:k+1,2),CPz(k:k+1,3),'--or');
  end
  i=i+1;
  CPz(:,1:3)=CP(i,j,:,1:3);
  plot3(CPz(k:k+1,1),CPz(k:k+1,2),CPz(k:k+1,3),'--or');
end
k=k+1;
for j = 1:length(CP(1,:,1,1))-1
  for i = 1:length(CP(:,1,1,1))-1
    plot3(CP(i:i+1,j,k,1),CP(i:i+1,j,k,2),CP(i:i+1,j,k,3),'--or');
    plot3(CP(i,j:j+1,k,1),CP(i,j:j+1,k,2),CP(i,j:j+1,k,3),'--or');
  end
  i=i+1;
  plot3(CP(i,j:j+1,k,1),CP(i,j:j+1,k,2),CP(i,j:j+1,k,3),'--or');
end
j=j+1;
for i = 1:length(CP(:,1,1,1))-1
  plot3(CP(i:i+1,j,k,1),CP(i:i+1,j,k,2),CP(i:i+1,j,k,3),'--or');
end