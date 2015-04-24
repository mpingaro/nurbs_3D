function i = findspan(u,U,n)
m=length(U);
% Note: due to rounding errors in u findspan can be ambiguous at knots.
% There's no problem for the inner knots but to get the last knot (special
% case) rounding range must be considered
eps=10e-10;
if (abs(u-U(n+1))<eps)  % special case: last knot, but works only for open
  i = n;                % knot vector!
  return
end
for i = 1:(m-1)
  if (u<U(i+1))
    return
  end
end  

error('u outside of U!')