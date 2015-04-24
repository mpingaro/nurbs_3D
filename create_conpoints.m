function create_conpoints(CP)
% plots control points without polygon

for k = 1:length(CP(1,1,:,1))
  for j = 1:length(CP(1,:,1,1))
    for i = 1:length(CP(:,1,1,1))
      plot3(CP(i,j,k,1),CP(i,j,k,2),CP(i,j,k,3),'or');
    end
  end
end