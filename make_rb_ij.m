function RB = make_rb_ij(rb,nu)
% reads a row vector of dofs and returns a matrix with the corresponding i
% and j indices and direction (1=x, 2=y)

for k = 1:length(rb)
  h=rb(k)/2;
  p=ceil(h);
  j=ceil(p/nu);
  i=p-(j-1)*nu;
  if (p~=h)  dir=1;
  else       dir=2;
  end
  
  RB(k,1)=i;
  RB(k,2)=j;
  RB(k,3)=dir;
end
