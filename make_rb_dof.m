function rb = make_rb_dof(RB,nu)
% reads a matrix with i and j indices and direction (1=x, 2=y) and returns
% the row vector of corresponding dofs

for k = 1:length(RB(:,1))
  i = RB(k,1);
  j = RB(k,2);
  cp = (j-1)*nu+i;
  if (RB(k,3)==1)
    rb(k) = 2*cp-1;
  else
    rb(k) = 2*cp;  
  end
end
